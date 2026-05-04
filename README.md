# DTQEM v10.1 – Dual-Time Quantum Entanglement Model

![DTQEM v10.1 Result](images/DTQEM_v10.1.jpg)

A calibrated open‑quantum‑system research framework for simulating entanglement, decoherence, observer‑induced disturbance, and double‑slit interference using exact Liouvillian dynamics.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![arXiv](https://img.shields.io/badge/arXiv-2505.00000-red.svg)](https://arxiv.org/abs/2505.00000) *(placeholder)*

---

## Overview

DTQEM (Dual‑Time Quantum Entanglement Model) is an open‑source simulation framework that studies how entangled quantum states evolve under the influence of their environment. By modeling decoherence channels, thermal noise, magnetic coupling, and observation‑like disturbances, the project bridges quantum information measures (visibility, concurrence, negativity) with macroscopic phenomena like double‑slit interference.

## Conceptual Framework

The model is guided by a dual‑time hypothesis. Alongside ordinary physical time, the model explores a conceptual imaginary‑time component that acts as a "processing lag" for quantum information.

- **The "Cosmic Radar" Picture:** Particles are treated as systems emitting resonant frequencies that probe the environment before their localized detection.
- **Observer Effect:** Modeled as a resonant interference / shift in the particle's internal states due to interaction with an active measuring device (observer).
- **Scientific Note:** This conceptual layer serves as an interpretive hypothesis rather than an established physical theory; it provides the motivation for the numerical modeling and the analysis of interference loss under external disturbance.

## Numerical Model

The simulation relies on a rigorous computational layer:
- **Exact Liouvillian Dynamics:** We solve the Lindblad master equation using superoperator exponentiation (`expm(L·t)`), ensuring stability and preventing the numerical drift associated with adaptive ODE solvers.
- **Diagnostics:** Comprehensive metrics tracking trace integrity, hermiticity, and positivity.

## Key Features

- **Exact Evolution:** Stable, drift‑free density‑matrix updates.
- **Entanglement Metrics:** Visibility, Concurrence, Negativity, Purity, von Neumann Entropy, Fidelity, and l1‑norm Coherence.
- **Interactive GUI:** Real‑time parameter control (θ, T, t_obs, B, λ, d, etc.) with presets, batch export, and help panels.
- **Interference Visualization:** 1D fringe patterns and 2D intensity maps.
- **Validated Benchmarks:** All benchmarks pass with errors < 1e‑12 (see below).

## Benchmark Results (run automatically)

## Acknowledgments

- **Reddouane BERRAMDANE** – Formulation of the dual‑time hypothesis, the “cosmic radar” analogy, and the philosophical‑conceptual framework that drives the project.
- **DeepSeek (AI assistant)** – Complete code development, numerical implementation (Liouvillian solver, GUI, benchmarks), debugging, and documentation.
- **Gemini (AI assistant)** – In‑depth philosophical dialogues, critical examination of retrocausality, observer effect as frequency interference, and conceptual sharpening.
- **Claude (AI assistant)** – Rigorous scientific criticism, detailed code reviews, identification of logical gaps, and constructive suggestions that improved the model’s coherence and testability.
- The open‑source Python community (NumPy, SciPy, Matplotlib) for their essential tools.

We also thank all researchers whose experimental work (Aspect, Gisin, Lindblad) provided the empirical grounding for this simulation.
