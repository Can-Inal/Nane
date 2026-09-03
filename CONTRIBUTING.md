# Contributing to Nane

Thank you for your interest in contributing to Nane! 🎉

Nane is still under active development, and feedback, bug reports, ideas, and
code contributions are all welcome.

## Questions and Discussions

If you have a question about using Nane, an API idea, or a general suggestion,
please use [GitHub Discussions](https://github.com/Can-Inal/Nane/discussions).

Before starting a new discussion, it is worth checking the
[documentation](https://can-inal.github.io/Nane/) and existing discussions in
case the topic has already come up.

## Issues

Please use [GitHub Issues](https://github.com/Can-Inal/Nane/issues) for bugs
and concrete feature requests.

For bug reports, useful information includes:

- what happened
- what you expected to happen
- a small reproducing example, if possible
- compiler and operating system
- Nane version or commit

For larger feature ideas, especially new numerical methods or changes to the
public API, opening an issue or discussion first is helpful. It gives me a
chance to think about the scope and design before implementation work begins.

Small fixes such as typos or obvious documentation corrections do not need a
separate issue first.

## Code Contributions

Pull requests are welcome.

If you are looking for something to work on, issues labeled `good first issue`
or `help wanted` are good places to start.

Please try to keep pull requests focused on one logical change and avoid
mixing unrelated refactoring, formatting, and new functionality.

Nane follows the formatting rules in `.clang-format` and generally uses:

- C++23
- the `nane` namespace
- `snake_case` names
- Allman braces
- small, composable interfaces
- Eigen types directly
- free functions for standalone numerical algorithms

For numerical contributions, I especially value clear correspondence between
the implementation and the underlying mathematics.

New functionality should normally include appropriate tests and public API
documentation.

## Numerical Methods

If you would like to contribute a numerical method, it is helpful to include
some context such as:

- the mathematical formulation
- expected accuracy or convergence behavior
- relevant references
- suitable test problems

This is not meant as a rigid checklist; it simply makes numerical changes
easier to review and understand.

## Commits

There is no strict commit-message convention.

Clear messages are appreciated, for example:

```text
Add classical Runge-Kutta method
Fix stage evaluation in implicit Runge-Kutta
Document fixed-point convergence criterion
```

I generally prefer imperative wording such as `Add`, `Fix`, `Remove`,
`Refactor`, or `Document`.

## Scope

Nane is intended to stay focused on numerical and scientific computing.

If you are unsure whether an idea belongs in Nane, feel free to start a
discussion before working on it.

## License

By contributing to Nane, you agree that your contributions will be licensed
under the terms of the project's [MIT License](LICENSE).
