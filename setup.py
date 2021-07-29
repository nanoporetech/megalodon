import sys
from setuptools import setup, Extension


if __name__ == "__main__":
    # Note that using setup_requires cython allows users to install megalodon
    # without first installing cython (as required when using cythonize)
    extra_compile_args = ["-std=c99"]
    if sys.platform == "darwin":
        extra_compile_args.append("-mmacosx-version-min=10.9")
        print("Using macOS clang args")
    ext_modules = [
        Extension(
            "megalodon.decode",
            sources=["megalodon/decode.pyx"],
            extra_compile_args=extra_compile_args,
            language="c",
        ),
    ]
    setup(
        use_pyscaffold=True,
        setup_requires=["setuptools>=38.3", "cython"],
        ext_modules=ext_modules,
    )
