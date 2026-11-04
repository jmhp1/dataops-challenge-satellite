"""Run the whole challenge in order:  python -m src.run_all

Each step only reads files written by earlier steps, so any step can also be
run on its own (python -m src.step3_maps), and an orchestrator (Prefect,
Dagster, Airflow) can wrap each step's main() as one task.
"""

import time

from src import step1_describe, step2_calibrate, step3_maps, step4_reproject, step5_merge, step6_weather

STEPS = [step1_describe, step2_calibrate, step3_maps, step4_reproject, step5_merge, step6_weather]

if __name__ == "__main__":
    for step in STEPS:
        start = time.time()
        print(f"\n=== {step.__name__.split('.')[-1]} ===")
        step.main()
        print(f"--- done in {time.time() - start:.0f} s")
