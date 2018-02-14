g++ -std=c++11 -stdlib=libc++ -v tests/test_graph.cpp graph.cpp tests/gtest_main.cpp -lgtest -lpthread -o tests/test_graph.out
./tests/test_graph.out
