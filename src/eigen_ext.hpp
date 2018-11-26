#ifndef pear_eigen_ext_hpp
#define pear_eigen_ext_hpp

#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/Core>

namespace pear {
template <typename T1, typename T2, typename T3>
T1 extract(const T2 &full, const T3 &ind) {
  int num_indices = ind.innerSize();
  T1 target(num_indices);
  for (int i = 0; i < num_indices; i++) {
    target[i] = full[ind[i]];
  }
  return target;
}

template <typename M> M load_csv(const std::string &path) {
  std::ifstream indata;
  indata.open(path);
  std::string line;
  std::vector<typename M::Scalar> values;
  uint64_t rows = 0;
  while (std::getline(indata, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ',')) {
      values.push_back(std::stod(cell));
    }
    ++rows;
  }
  return Eigen::Map<M>(values.data(), rows, values.size() / rows);
}

template <typename M> void write_csv(const std::string &path, M &matrix) {
  std::ofstream file(path.c_str());

  for (int i = 0; i < matrix.rows(); i++) {
    for (int j = 0; j < matrix.cols(); j++) {
      std::string str = std::to_string(matrix(i, j));
      if (j + 1 == matrix.cols()) {
        file << str;
      } else {
        file << str << ',';
      }
    }
    file << '\n';
  }
}

} // namespace pear

#endif
