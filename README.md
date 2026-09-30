# Schwinger–Dyson Equation Generator for Multi-Matrix Models

Java programs that generate the large-`N` Schwinger–Dyson equations of multi-matrix models, reduce them by symmetry, and export them as runnable Mathematica code. Also programs that find all multi-matrix potentials with sufficient symmetry. Written for:

> M. Khalkhali, N. Pagliaroli, A. Parfeni, B. Smith, *Bootstrapping the Critical Behavior of Multi-Matrix Models*, **J. High Energ. Phys. 2025, 158** — [arXiv:2409.07565](https://arxiv.org/abs/2409.07565)

## How it works

1. **Generate symmetric potentials** (batch mode): build every model from words up to a given length and keep those unchanged under any relabelling of the matrices.
2. **Enumerate** every word (moment) up to a truncation length and derive its Schwinger–Dyson equation.
3. **Reduce** each moment to a canonical representative under cyclic shifts, transposition, matrix relabelling, and optionally sign flips.
4. **Export** the reduced system to Mathematica with the `Solve` call included.

## Generating symmetric potentials

A potential is a set of trace words (`ABAB` = tr(ABAB)), all with coupling `g`. Candidate terms are all words up to `maxLengthWords` in which each matrix appears an even number of times (any count in `*Odd`), taken up to cyclic shift. Every combination of up to `maxNumTerms` terms is tested, and a model is kept only if relabelling the matrices maps its set of terms to itself. Single-matrix quadratic and quartic terms (cubic in `*Odd`) are added automatically and excluded from the search.

Example, two matrices: `{AABB, ABAB}` passes; `{AAAABB}` fails, since `A↔B` gives `AABBBB`; `{AAAABB, AABBBB}` passes.

Feasible sizes (matrices, word length, terms): below (4, 6, 4), (3, 6, 6), (2, 8, 8).

## Files

| File | Purpose |
|---|---|
| `GeneralSDEMathematica.java` | **Single model.** Takes one potential you specify in `extraWords`, generates its Schwinger–Dyson system to the chosen truncation, and writes Mathematica code that solves it and lists the moments left undetermined (the ones the bootstrap must bound). Supports sign-flip symmetry. |
| `GeneralSDEMathematicaFull.java` | **Batch mode.** Generates every symmetric potential within the size limits, then writes one Schwinger–Dyson system per potential as a Mathematica list, with a helper that solves each and extracts the same moment, so a whole family of models is scanned in one evaluation. |
| `GeneralSDEMathematicaOdd.java` | Single model, allowing words where a matrix appears an odd number of times, for potentials with odd-degree terms (e.g. cubic). |
| `GeneralSDEMathematicaFullOdd.java` | Batch mode for odd-degree potentials; the automatic single-matrix term is cubic instead of quartic. |
| `symmetricModels.java` | Only the potential search: writes every symmetric potential to `output2.txt` without generating equations. Use it to size the search space before a batch run. |
| `symmetricModelsOdd.java` | Potential search for odd-degree potentials. |

## Usage

Requires a JDK and Mathematica. Set parameters at the top of `main`, then:

```bash
javac GeneralSDEMathematica.java && java GeneralSDEMathematica
```

Paste `output.txt` into Mathematica and evaluate. Moments appear as `x1, x2, …` (mapping to words printed while running); couplings `t1, t2, …` are set to `g`.

| Parameter | Meaning |
|---|---|
| `lNum` | Number of matrices. |
| `extraWords` | Single model: the potential, e.g. `{"ABAB", "AABB"}`. Batch: its longest word sets the term length. |
| `maxLengthofWordsGeneratingSDEs` | Truncation length. ≤5 instant, ≤7 under a minute, ≥9 can be very slow. |
| `maxNumTerms` / `maxLengthWords` | Batch and potential search: max terms per model, max term length. |
| `fullSymmetry` / `flipSymmetry` | Relabelling group (full symmetric or cyclic only); sign-flip symmetry. |

## Known issues

- The `GeneralSDE*` files largely duplicate each other; parameters are edited in source.

