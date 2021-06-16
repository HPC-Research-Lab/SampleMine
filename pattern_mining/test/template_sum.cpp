#include <iostream>
using namespace std;

template<int...>
struct sum;

template<int s>
struct sum<s> {
  enum {value = s};
};

template<int s, int... others>
struct sum<s, others...> {
  enum {value = s + sum<others...>::value };
};

template<int s1, int s2, int... others>
struct sum<s1, s2, others...> {
  enum {value = s1 + sum<s2, others...>::value };
};

int main() {
  cout << sum<1,2,3,4,5>::value << endl;
  return 0;

}