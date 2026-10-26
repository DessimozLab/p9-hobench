"""HOBench tool wrapper (DESIGN §8).

* ``python -m wrapper.run``     — run OrthoFinder / FastOMA (or derive the of_flat baseline) on a
  panel; work dirs on local NVMe, outputs + manifest.json in runs/<panel>/<mode>/<tool>_<variant>/.
* ``python -m wrapper.subset``  — family-neighbourhood subset panels (SYNTHESIS D5).
* ``python -m wrapper.perturb`` — robustness perturbations of a panel (SYNTHESIS E5).
* ``scripts/run_full_YK.sh``    — queue of the full-mode Y/K runs.

Blind to reference families: no module reads families/.
"""
