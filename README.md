<p align="center">
    <img src="docs/assets/nane_logo.png" alt="Nane" width="480">
</p>

<p align="center">
    <a href="https://github.com/Can-Inal/Nane/actions/workflows/build.yml">
        <img src="https://img.shields.io/github/actions/workflow/status/Can-Inal/Nane/build.yml?branch=main&style=for-the-badge&logo=github&logoColor=white&label=Build" alt="Build">
    </a>
    <a href="https://can-inal.github.io/Nane/">
        <img src="https://img.shields.io/badge/Documentation-Doxygen-2C4AA8?style=for-the-badge&logo=doxygen&logoColor=white" alt="Documentation">
    </a>
    <img src="https://img.shields.io/badge/C%2B%2B-23-00599C?style=for-the-badge&logo=cplusplus&logoColor=white" alt="C++23">
    <img src="https://img.shields.io/badge/CMake-3.20%2B-064F8C?style=for-the-badge&logo=cmake&logoColor=white" alt="CMake">
    <a href="LICENSE">
        <img src="https://img.shields.io/github/license/Can-Inal/Nane?style=for-the-badge" alt="License">
    </a>
</p>

Nane is a modern C++23 scientific computing library for building and exploring numerical algorithms with a clean and direct interface.

It is designed to keep numerical code understandable, lightweight, and easy to extend.

> **Nane is currently under active development.**

## Build

```bash
cmake -S . -B build -DNANE_BUILD_EXAMPLES=ON
cmake --build build --parallel
```

Examples are available in the [`examples`](examples/) directory.

## Documentation

The API documentation is generated with Doxygen and published through GitHub Pages.

[Open the Nane documentation](https://can-inal.github.io/Nane/)

## Contributing

Contributions, bug reports, feature suggestions, and discussions are welcome.

See [CONTRIBUTING.md](CONTRIBUTING.md) for more information.

## License

Nane is available under the [MIT License](LICENSE).
