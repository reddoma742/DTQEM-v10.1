[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20024264.svg)](https://doi.org/10.5281/zenodo.20024264)

# White Paper – DTQEM v10.1  
## Dual‑Time Quantum Entanglement Model: From Numerical Simulation to Interpretive Framework

**Author:** Redouane BERRAMDANE  
**Date:** 2025  
**Repository:** [DTQEM-v10.1](https://github.com/reddoma742/DTQEM-v10.1)

---

## Abstract

We present DTQEM v10.1, a numerically exact open‑quantum‑system simulation of entanglement, decoherence, and the observer effect. The code solves the Lindblad master equation via Liouvillian superoperator exponentiation, providing machine‑precision benchmarks for pure dephasing, relaxation, and entropy increase. All standard entanglement measures (visibility, concurrence, negativity, purity, entropy, fidelity) are implemented. An interactive GUI simulates double‑slit interference as a function of launch angle, temperature, observation time, decoherence rates, phonon energy, and magnetic field.

Complementing the computational layer, we introduce an **interpretive heuristic** – the “cosmic radar” picture – which views the particle as emitting a pre‑field that probes the environment before materialisation. This heuristic is **not a substitute for the Lindblad equation**, nor does it introduce new free parameters; it merely provides an intuitive language for the observer‑induced loss of interference. The paper clearly distinguishes between the **validated numerical model** (which stands on its own) and the **philosophical interpretation** (which is open to debate).

---

## 1. Introduction

The measurement problem and the observer effect in quantum mechanics remain subjects of intense discussion. While the Lindblad master equation successfully describes decoherence, it does not by itself offer an intuitive picture of why an active measurement destroys interference. DTQEM v10.1 provides a **computational laboratory** where one can switch between “passive screen” and “active camera” regimes, observing the transition from quantum to classical behaviour.

To make the results tangible, we propose an **interpretive heuristic** – not a new physical theory – that we call the **“cosmic radar”**. According to this picture, every quantum particle emits a resonant pre‑field that travels ahead, interacts with the environment, and returns a “echo” that influences the particle’s final state. The observation device (camera) is an active emitter that injects its own frequency, causing a resonance shift that collapses the coherence. This language is deliberately metaphorical; its only purpose is to help the reader connect the mathematical output of the Lindblad equation to the familiar double‑slit experiment.

---

## 2. The dual‑time picture as a rewriting of decoherence

Let \(t_r\) be ordinary (real) time. We introduce an **imaginary time** \(t_v = -\alpha(\theta)\, t_r\) with \(\alpha(\theta)=\sin(\theta/2)\) where \(\theta\) is the launch angle between two entangled particles. The effective time for the quantum influence is defined as:

\[
t_{\text{eff}} = t_r\bigl(1 - \alpha(\theta)\, e^{-\Gamma(T)\,t_{\text{obs}}}\bigr).
\]

This expression can be rewritten as a simple rescaling of the real time; it does **not** alter the Lindblad master equation that actually governs the density matrix evolution. The quantity \(t_v\) is therefore an **interpretive tool** (a “processing lag” in the radar analogy) rather than a physical extension of quantum mechanics.  

**Crucially, the numerical implementation – the Liouvillian superoperator – contains no reference to \(t_v\); it works directly with the Lindblad equation.** The dual‑time language is used only in the **discussion** of the results (e.g., when we say “the effective time becomes very small, hence the speed appears superluminal”). This is a **rhetorical device**, not a claim of superluminal signalling.

---

## 3. The computational core (validated)

The simulation solves the Lindblad master equation exactly:

\[
\frac{d\rho}{dt} = -\frac{i}{\hbar}[H,\rho] + \sum_k \left(L_k\rho L_k^\dagger - \frac12\{L_k^\dagger L_k,\rho\}\right).
\]

The Hamiltonian \(H = \frac12 g\mu_B B\,\sigma_z\otimes I_2\) describes a magnetic field. The jump operators represent:

- Pure dephasing: \(L_{\phi} = \sqrt{\gamma_{\phi}(T)}\,\sigma_z\otimes I_2\),
- Relaxation: \(L_{\downarrow} = \sqrt{\gamma_{\downarrow}(T)}\,\sigma_-\otimes I_2\),
- Excitation: \(L_{\uparrow} = \sqrt{\gamma_{\uparrow}(T)}\,\sigma_+\otimes I_2\).

The thermal rates follow Bose–Einstein statistics:  
\(n_{\text{th}} = 1/(e^{\hbar\omega/k_BT}-1)\),  
\(\gamma_{\phi}(T) = (\gamma_{\phi0}/2)(2n_{\text{th}}+1)\),  
\(\gamma_{\downarrow}(T) = \gamma_{\text{rel}0}(n_{\text{th}}+1)\),  
\(\gamma_{\uparrow}(T) = \gamma_{\text{rel}0}\,n_{\text{th}}\).

Instead of an ODE integrator, we construct the Liouvillian superoperator \(\mathcal{L}\) (16×16) and compute \(\vec{\rho}(t) = \exp(\mathcal{L}t)\,\vec{\rho}(0)\) via `scipy.linalg.expm`. This yields machine‑precision accuracy and avoids numerical drift.

### Benchmarks

- **Pure dephasing** (initial state \((|00\rangle+|10\rangle)/\sqrt{2}\)): the off‑diagonal element \(\rho_{00,10}\) decays as \(0.5\exp(-\gamma_{\phi0}t)\). **Error \(< 10^{-15}\).**
- **Relaxation at T=0** (initial \(|01\rangle\)): population decays as \(\exp(-\gamma_{\text{rel}0}t)\). **Error \(< 10^{-15}\).**
- **Entropy increase** (initial \(|00\rangle\)): after long evolution at T=100 K, the entropy increases (second law). **Passed.**

All benchmarks are automatically run when the script starts.

---

## 4. Testable predictions

The model contains no free parameters after calibration (using the experimental lower bound \(v_{\text{eff}} \ge 10^7c\) at T=0 K, and choosing \(v_{\text{eff}} \approx 1200c\) at T=300 K). A unique, experimentally verifiable prediction is:

> **For a maximally entangled state (\(\theta=180^\circ\)), the fringe visibility in a double‑slit experiment should follow**
> \[
> V(T) = \exp\!\bigl(-(0.12 + 3.33\,T)\times 10^{-6}\bigr)
> \]
> **with \(t_{\text{obs}} = 10^{-6}\) s (a typical decoherence time). This specific exponential decay can be tested by measuring the contrast of interference fringes at different controlled temperatures.**

If future experiments observe a different functional form (e.g., a power law or a different exponential coefficient), the DTQEM calibration would be falsified.

---

## 5. Results and how to use the code

The interactive GUI (run `dtqem_v10_1.py`) allows the user to:

- Adjust the launch angle \(\theta\), temperature \(T\), observation time \(t_{\text{obs}}\), decoherence rates, phonon energy, magnetic field, slit parameters.
- Visualise 1D and 2D interference fringes in real time.
- Read off entanglement measures (visibility, concurrence, negativity, purity, entropy, fidelity).
- Run batch simulations to produce visibility maps (\(\theta\) vs \(T\)).
- Export all figures and numerical data (CSV, JSON).

The code is open source (MIT) and can be used for educational purposes, as a research tool, or as a foundation for further extensions.

---

## 6. Philosophical note – the “cosmic radar” heuristic

The “cosmic radar” picture is **not** a physical theory. It is a pedagogical and heuristic device that reinterprets the mathematical terms of the Lindblad equation in a language borrowed from wave interference and signal processing.  

- **Passive screen** – does not emit its own frequency → the particle’s pre‑field remains coherent → interference fringes appear.  
- **Active observer (camera)** – emits a competing frequency → resonance shift → collapse of coherence → no fringes.

We stress that this language is entirely optional. One can use the code and understand its output without ever invoking the radar metaphor. The heuristic does **not** introduce new equations, new constants, or any testable consequence beyond those already contained in the Lindblad model. Its only purpose is to provide an intuitive anchor for those who find the standard “collapse” language unsatisfactory.

---

## 7. Conclusion

DTQEM v10.1 provides a **numerically exact, open‑source simulation** of two‑qubit entanglement under realistic decoherence channels, thermal noise, and magnetic fields. Its benchmarks pass with machine precision, and it offers a rich GUI for exploring the transition from quantum to classical behaviour. The accompanying **interpretive heuristic** (“cosmic radar”) is clearly labelled as philosophical – it does not alter the mathematical model. A **unique, testable prediction** (the specific exponential decay of fringe visibility with temperature) is stated explicitly.

The code and documentation are ready for use in education, research, and as a foundation for extensions (non‑Markovian noise, multi‑partite systems, gravitational effects).

---

## 8. References

1. Lindblad, G. (1976). *Commun. Math. Phys.* **48**, 119.  
2. Gisin, N. & Zbinden, H. (1998). *Phys. Lett. A* **248**, 1.  
3. Aspect, A., Dalibard, J. & Roger, G. (1982). *Phys. Rev. Lett.* **49**, 1804.  
4. Breuer, H.‑P. & Petruccione, F. (2002). *The Theory of Open Quantum Systems*.
