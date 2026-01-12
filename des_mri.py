

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd

WORK_START_H_DEFAULT = 8
WORK_END_H_DEFAULT = 17

@dataclass
class Inputs:
    lambda_per_day: float
    mu_duration_min: float
    sigma_duration_min: float
    type2_daily_arrivals: np.ndarray          # integer counts per observed day
    type2_durations_min: np.ndarray           # empirical durations in minutes
    slot_minutes_type1: int
    slot_minutes_type2: int
    work_minutes: int = 540                   # 9 hours * 60
    seed: int = 123

def _is_workday(day_index: int) -> bool:
    """Assume day_index=0 is Monday. Workdays are Mon-Fri."""
    dow = day_index % 7  # 0..6
    return dow <= 4      # 0..4 => Mon..Fri

def _next_workday(day_index: int) -> int:
    d = day_index + 1
    while not _is_workday(d):
        d += 1
    return d

def _generate_calls_for_day(inputs: Inputs, rng: np.random.Generator, day_index: int) -> pd.DataFrame:
    """
    Generate call times within the working day for both types.
    Returns DataFrame with columns: day_index, call_time_min, patient_type
    """
    W = inputs.work_minutes

    # Type 1: Poisson process over working minutes
    rate_per_min = inputs.lambda_per_day / W
    t = 0.0
    t1 = []
    while True:
        t += rng.exponential(1.0 / rate_per_min)
        if t > W:
            break
        t1.append(t)

    # Type 2: sample daily count from empirical and place uniformly
    n2 = int(rng.choice(inputs.type2_daily_arrivals))
    t2 = np.sort(rng.uniform(0, W, size=n2)) if n2 > 0 else np.array([])

    calls = []
    for x in t1:
        calls.append((day_index, float(x), "Type 1"))
    for x in t2:
        calls.append((day_index, float(x), "Type 2"))

    if not calls:
        return pd.DataFrame(columns=["day_index", "call_time_min", "patient_type"])

    df = pd.DataFrame(calls, columns=["day_index", "call_time_min", "patient_type"])
    df = df.sort_values(["call_time_min", "patient_type"], kind="mergesort").reset_index(drop=True)
    return df

def _candidate_slot_starts(work_minutes: int, slot_len: int) -> np.ndarray:
    """Slot starts (minutes since 08:00) such that slot fits within day."""
    if slot_len <= 0 or slot_len > work_minutes:
        return np.array([], dtype=int)
    # start times at multiples of slot_len, from 0 to W-slot_len inclusive
    return np.arange(0, work_minutes - slot_len + 1, slot_len, dtype=int)

def _book_earliest_slot(
    schedule: Dict[Tuple[int, int], set],
    machines: List[int],
    earliest_day: int,
    patient_type: str,
    inputs: Inputs,
    policy: str
) -> Tuple[int, int, int]:
    """
    Book earliest available slot for the patient under given policy.

    schedule key: (machine_id, day_index) -> set(start_times)
    Returns (appointment_day, machine_id, slot_start_min)
    """
    slot_len = inputs.slot_minutes_type1 if patient_type == "Type 1" else inputs.slot_minutes_type2
    starts = _candidate_slot_starts(inputs.work_minutes, slot_len)
    if len(starts) == 0:
        raise ValueError("No feasible slots (slot length too long for workday).")

    d = earliest_day
    while True:
        if not _is_workday(d):
            d = _next_workday(d - 1)  # jump to next workday
            continue

        # determine eligible machines
        if policy == "old":
            eligible = [0] if patient_type == "Type 1" else [1]
        elif policy == "new":
            eligible = machines
        else:
            raise ValueError("policy must be 'old' or 'new'")

        # Scan slots in chronological order across eligible machines.
        # We select the earliest (day, start_time), breaking ties by machine_id.
        best = None
        for m in eligible:
            booked = schedule.setdefault((m, d), set())
            for s in starts:
                if s not in booked:
                    cand = (d, s, m)
                    if best is None or cand < best:  # tuple order: day, start, machine
                        best = cand
                    break  # earliest free slot on this machine for day d

        if best is not None:
            appt_day, start_min, machine_id = best
            schedule[(machine_id, appt_day)].add(start_min)
            return appt_day, machine_id, int(start_min)

        # no slot on this day, go to next workday
        d = _next_workday(d)

def _simulate_day_operations(
    assigned_today: pd.DataFrame,
    inputs: Inputs,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """
    Given all patients scheduled for a specific day with assigned machine + slot start,
    simulate actual scan durations and compute per-machine busy time and overtime.
    Returns per-machine metrics for that day.
    """
    if assigned_today.empty:
        # return 2 machines with zeros
        return pd.DataFrame({
            "machine_id": [0, 1],
            "n_scans": [0, 0],
            "busy_min": [0.0, 0.0],
            "finish_min": [0.0, 0.0],
            "overtime_min": [0.0, 0.0],
            "utilization": [0.0, 0.0],
        })

    df = assigned_today.copy()

    # sample actual durations
    dur = np.zeros(len(df), dtype=float)
    is_t1 = (df["patient_type"].values == "Type 1")
    n1 = int(is_t1.sum())
    n2 = len(df) - n1
    if n1 > 0:
        d1 = rng.normal(inputs.mu_duration_min, inputs.sigma_duration_min, size=n1)
        d1 = np.where(d1 <= 0, 0.1, d1)
        dur[is_t1] = d1
    if n2 > 0:
        dur[~is_t1] = rng.choice(inputs.type2_durations_min, size=n2, replace=True)

    df["actual_duration_min"] = dur

    out_rows = []
    for m in [0, 1]:
        dm = df[df["machine_id"] == m].sort_values("slot_start_min")
        if dm.empty:
            out_rows.append((m, 0, 0.0, 0.0, 0.0, 0.0))
            continue

        # Machine processes sequentially; if a scan runs long, it pushes the next start (spillover).
        t = 0.0
        busy = 0.0
        for _, r in dm.iterrows():
            start = float(r["slot_start_min"])
            # actual start is max(planned slot start, previous completion)
            t = max(t, start)
            t = t + float(r["actual_duration_min"])
            busy += float(r["actual_duration_min"])

        finish = t
        overtime = max(0.0, finish - inputs.work_minutes)
        util = busy / inputs.work_minutes
        out_rows.append((m, int(len(dm)), busy, finish, overtime, util))

    daily = pd.DataFrame(out_rows, columns=["machine_id", "n_scans", "busy_min", "finish_min", "overtime_min", "utilization"])
    return daily

def run_simulation(
    inputs: Inputs,
    n_workdays: int = 250,
    policy: str = "old",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run simulation for n_workdays workdays (not calendar days).
    We simulate a calendar where day_index=0 is Monday, and we continue until we have simulated
    n_workdays workdays worth of operations.

    Returns:
      patients_df: call + assignment + realized durations
      daily_machine_df: per workday, per machine metrics
    """
    rng = np.random.default_rng(inputs.seed)
    machines = [0, 1]
    schedule: Dict[Tuple[int, int], set] = {}  # (machine, day) -> set(slot_starts)
    patients_records = []
    daily_records = []

    workdays_done = 0
    day_index = 0

    while workdays_done < n_workdays:
        if not _is_workday(day_index):
            day_index += 1
            continue

        # Generate calls for this day and book appointments
        calls = _generate_calls_for_day(inputs, rng, day_index)

        for _, c in calls.iterrows():
            call_day = int(c["day_index"])
            call_time = float(c["call_time_min"])
            ptype = str(c["patient_type"])

            earliest = _next_workday(call_day)
            appt_day, machine_id, slot_start = _book_earliest_slot(
                schedule=schedule,
                machines=machines,
                earliest_day=earliest,
                patient_type=ptype,
                inputs=inputs,
                policy=policy
            )

            slot_len = inputs.slot_minutes_type1 if ptype == "Type 1" else inputs.slot_minutes_type2

            patients_records.append({
                "call_day": call_day,
                "call_time_min": call_time,
                "patient_type": ptype,
                "appointment_day": appt_day,
                "machine_id": machine_id,
                "slot_start_min": slot_start,
                "slot_len_min": int(slot_len),
            })

        # Simulate operations for THIS day (patients whose appointment_day == day_index)
        patients_df_partial = pd.DataFrame(patients_records)
        assigned_today = patients_df_partial[patients_df_partial["appointment_day"] == day_index].copy()

        daily_m = _simulate_day_operations(assigned_today, inputs, rng)
        daily_m["day_index"] = day_index
        daily_records.append(daily_m)

        workdays_done += 1
        day_index += 1

    patients_df = pd.DataFrame(patients_records)
    # waiting time in working days (count of workdays between call_day and appointment_day)
    # Because we simulate a week calendar, compute it accurately:
    def working_days_between(a: int, b: int) -> int:
        if b <= a:
            return 0
        cnt = 0
        for d in range(a+1, b+1):
            if _is_workday(d):
                cnt += 1
        return cnt

    if not patients_df.empty:
        patients_df["wait_workdays"] = [
            working_days_between(int(cd), int(ad))
            for cd, ad in zip(patients_df["call_day"].values, patients_df["appointment_day"].values)
        ]

    daily_machine_df = pd.concat(daily_records, ignore_index=True)

    # Attach actual durations back to patients_df for those simulated days.
    # We re-simulated operations day-by-day; to keep the final patient-level durations,
    # we can regenerate them consistently by re-running with the same seed and capturing per-day assignments.
    # For KPI purposes (wait times/overtime/utilization), patient-level durations are optional.
    # If needed, you can extend the simulation to store actual_duration_min inside _simulate_day_operations.

    return patients_df, daily_machine_df

def load_inputs_from_csv(
    type1_params_csv: str,
    type2_daily_arrivals_csv: str,
    type2_durations_csv: str,
    metadata_csv: str,
    seed: int = 123
) -> Inputs:
    p1 = pd.read_csv(type1_params_csv)
    meta = pd.read_csv(metadata_csv)

    def get_param(name: str) -> float:
        row = p1[p1["parameter"] == name]
        if row.empty:
            raise ValueError(f"Missing parameter '{name}' in {type1_params_csv}")
        return float(row["value"].iloc[0])

    def get_meta(key: str) -> str:
        row = meta[meta["key"] == key]
        if row.empty:
            raise ValueError(f"Missing key '{key}' in {metadata_csv}")
        return str(row["value"].iloc[0])

    lambda_per_day = get_param("lambda_per_day")
    mu_duration_min = get_param("mu_duration_min")
    sigma_duration_min = get_param("sigma_duration_min")

    slot1 = int(float(get_meta("slot_minutes_type1")))
    slot2 = int(float(get_meta("slot_minutes_type2")))
    work_minutes = int(float(get_meta("work_minutes"))) if "work_minutes" in set(meta["key"]) else 540

    d2_daily = pd.read_csv(type2_daily_arrivals_csv)["daily_arrivals"].astype(int).to_numpy()
    d2_dur = pd.read_csv(type2_durations_csv).iloc[:, 0].astype(float).to_numpy()

    return Inputs(
        lambda_per_day=lambda_per_day,
        mu_duration_min=mu_duration_min,
        sigma_duration_min=sigma_duration_min,
        type2_daily_arrivals=d2_daily,
        type2_durations_min=d2_dur,
        slot_minutes_type1=slot1,
        slot_minutes_type2=slot2,
        work_minutes=work_minutes,
        seed=seed
    )
