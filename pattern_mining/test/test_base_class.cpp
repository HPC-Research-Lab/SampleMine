#include <iostream> 

class Base {
public:
  virtual void work() = 0;
};

class D1 : public Base {
public:
  void work() {
    std::cout << "d1" << std::endl;
  }
};

void func(Base& t) {
  t.work();
}

int main() {
  auto t = D1();
  func(t);
  return 0;
}