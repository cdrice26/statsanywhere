cmake -S . -B build
cmake --build build --target test_runner --config Release
cd build && ctest -V
cd ..
