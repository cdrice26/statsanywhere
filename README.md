# StatsAnywhere
StatsAnywhere is an experimental, incomplete statistics & scientific computing library implemented in C. Currently, it supports some common linear algebra computations, special functions, and probability distribution functions. Special functions are computed using Gaussian quadrature.

## Why Do This When So Many Scientific Libraries Already Exist?
The intention was mostly for this to be a learning project into the underpinnings of probability distributions as well as low-level coding. It also refreshed some of my linear algebra and calculus knowledge along the way. It was honestly a fun challenge - while I did ask AI to get me unstuck a few times, particularly in how to calculate the values for the Jacobi matrix in the Golub-Welsch algorithm and in converting from a general QR decomposition process to a tridiagonal one, by the end of it I was doing most of the work independently and just using it to check my work for silly mistakes (and set up CMake scaffolding), so it's safe to say I learned a fair bit from this project! The only actual practical reason outside personal learning to create this is that implmenting it in C allows it to be imported by virtually any programming language, adding scientific and statistical capabilities to languages that currently lack robust libraries for these tasks.

## Should I Use It?
Probably not. It's not nearly as battle-tested as most scientific libraries, and some functions produce noticeably less precise output. Additionally, since the goal was just to get probability functions working, it's not really a full scientific library anyway. Feel free to download it and play around with it, but it's definitely not production-ready.

## Is It Documented?
Yes, the library has documentation comments and a complete API documentation is available [here](./API.md).

## How Do I Build It?
StatsAnywhere uses CMake for building, which can be complicated. For a simple way to build the app, use the `scripts/build.py` script, then link against the produced static library like any other C library.

## Running Tests
Testing uses [Unity](https://www.throwtheswitch.org/unity). Before running tests, it is necessary to install unity into the project. To do so (in a POSIX-compliant shell) run the following:
```bash
mkdir -p tests/unity/src tests/unity/include
curl -L https://raw.githubusercontent.com/ThrowTheSwitch/Unity/master/src/unity.c -o tests/unity/src/unity.c
curl -L https://raw.githubusercontent.com/ThrowTheSwitch/Unity/master/src/unity.h -o tests/unity/include/unity.h
curl -L https://raw.githubusercontent.com/ThrowTheSwitch/Unity/master/src/unity_internals.h -o tests/unity/include/unity_internals.h
```
You can then run tests with `scripts/test.py`. Note that if you change any tests and need to regenerate the reference file, you'll need to edit `scripts/generate_reference.py` to produce the correct file, then create a `venv`, install scipy, and run `scripts/generate_reference.py` inside it.
