#include <vector>
#include <set>
#include <iostream>

int main() {


  std::vector<std::vector<int>> t;

  for (int i = 0; i < 4; i++) {
    t.emplace_back(std::vector<int>(1024*1024*1024));
  }

  std::cout << "before clear" << std::endl;
  std::cin.get();
  t.clear();
  std::cout << "called clear" << std::endl;
  std::cin.get();



}