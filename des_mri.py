

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
    type2_daily_arrivals: np.ndarray
    type2_durations_min: np.ndarray
    slot_minutes_type1: int
    slot_minutes_type2: int
    work_minutes: int = 540
    seed: int = 123

def _is_workday(day_index: int) -> bool:

    dow = day_index % 7
    return dow <= 4

def _next_workday(day_index: int) -> int:
    d = day_index + 1
    while not _is_workday(d):
        d += 1
    return d

def _generate_calls_for_day(inputs: Inputs, rng: np.random.Generator, day_index: int) -> pd.DataFrame:

    W = inputs.work_minutes

    rate_per_min = inputs.lambda_per_day / W
    t = 0.0
    t1 = []
    while True:
        t += rng.exponential(1.0 / rate_per_min)
        if t > W:
            break
        t1.append(t)

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
    if slot_len <= 0 or slot_len > work_minutes:
        return np.array([], dtype=int)
    return np.arange(0, work_minutes - slot_len + 1, slot_len, dtype=int)

def _book_earliest_slot(
    schedule: Dict[Tuple[int, int], set],
    machines: List[int],
    earliest_day: int,
    patient_type: str,
    inputs: Inputs,
    policy: str
) -> Tuple[int, int, int]:

    slot_len = inputs.slot_minutes_type1 if patient_type == "Type 1" else inputs.slot_minutes_type2
    starts = _candidate_slot_starts(inputs.work_minutes, slot_len)
    if len(starts) == 0:
        raise ValueError("No feasible slots (slot length too long for workday).")

    d = earliest_day
    while True:
        if not _is_workday(d):
            d = _next_workday(d - 1)
            continue

        if policy == "old":
            eligible = [0] if patient_type == "Type 1" else [1]
        elif policy == "new":
            eligible = machines
        else:
            raise ValueError("policy must be 'old' or 'new'")


        best = None
        for m in eligible:
            booked = schedule.setdefault((m, d), set())
            for s in starts:
                if s not in booked:
                    cand = (d, s, m)
                    if best is None or cand < best:
                        best = cand
                    break

        if best is not None:
            appt_day, start_min, machine_id = best
            schedule[(machine_id, appt_day)].add(start_min)
            return appt_day, machine_id, int(start_min)

        d = _next_workday(d)

def _simulate_day_operations(
    assigned_today: pd.DataFrame,
    inputs: Inputs,
    rng: np.random.Generator,
) -> pd.DataFrame:

    if assigned_today.empty:
        return pd.DataFrame({
            "machine_id": [0, 1],
            "n_scans": [0, 0],
            "busy_min": [0.0, 0.0],
            "finish_min": [0.0, 0.0],
            "overtime_min": [0.0, 0.0],
            "utilization": [0.0, 0.0],
        })

    df = assigned_today.copy()

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

        t = 0.0
        busy = 0.0
        for _, r in dm.iterrows():
            start = float(r["slot_start_min"])
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

    rng = np.random.default_rng(inputs.seed)
    machines = [0, 1]
    schedule: Dict[Tuple[int, int], set] = {}
    patients_records = []
    daily_records = []

    workdays_done = 0
    day_index = 0

    while workdays_done < n_workdays:
        if not _is_workday(day_index):
            day_index += 1
            continue

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

        patients_df_partial = pd.DataFrame(patients_records)
        assigned_today = patients_df_partial[patients_df_partial["appointment_day"] == day_index].copy()

        daily_m = _simulate_day_operations(assigned_today, inputs, rng)
        daily_m["day_index"] = day_index
        daily_records.append(daily_m)

        workdays_done += 1
        day_index += 1

    patients_df = pd.DataFrame(patients_records)

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
