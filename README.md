# Discrete Variable Representation (DVR) Method

The DVR method discretizes continuous parameters to numerically solve systems, with basis functions localized around discrete grid points. Operators like the coordinate are diagonal in this representation, enabling efficient transformations between finite basis sets and finite grids spanning the system of interest.

### Key Features
- **Localized Basis Functions**: Approximate coordinate operators are diagonal at DVR points.
- **Efficient Representation**: Transforms between discrete basis sets and finite spatial grids.
- **Applications**: Numerical solutions of wavefunctions for bound potentials, e.g., harmonic oscillators.

### Example
For a harmonic oscillator, the DVR method represents continuous basis \( \delta(x) \) of Hamiltonian eigenstates in an approximate grid basis \( \theta(i) \), allowing efficient numerical solutions.
