#include "db.h"

#include <fstream>
#include <sstream>
#include <typeinfo>

#include "util.h"

namespace euler::db {

  static int db_path_num = -1;

  template <class key_type, class value_type>
  int MyKV<key_type, value_type>::db_id = 0;

  template <class key_type, class value_type>
  MyKV<key_type, value_type>::MyKV(int nc) : ncols(nc), data_size(0) {
    DBPath = std::string(std::getenv("DB_PATH"));
    assert(DBPath != "");
    if (db_path_num == -1) {
      srand(time(NULL));
      db_path_num = rand();
      DBPath += "/" + std::to_string(db_path_num);
      std::cout << "db path: " << DBPath << std::endl;
      assert(mkdir(DBPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0);
    }
    else {
      DBPath += "/" + std::to_string(db_path_num);
    }

    DBPath += std::string("/MyKV_path_") + typeid(key_type).name() +
      std::to_string(db_id++);

    struct stat sb;

    if (stat(DBPath.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
      int s = system(("rm -r " + DBPath).c_str());
    }

    assert(mkdir(DBPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0);
    nfiles = 0;
  }

  template <>
  void MyKV<std::string>::put(const void* a, size_t len) {
    keys[std::to_string(nfiles)] = nfiles;
    int fd = open((DBPath + "/" + std::to_string(nfiles++) + ".dat").c_str(),
      O_WRONLY | O_CREAT, (mode_t)0600);
    if (fd == -1) {
      perror("store data failed");
      exit(EXIT_FAILURE);
    }
    ssize_t r = write(fd, a, len);
    if (r < 0) {
      perror("db error 1");
      exit(-1);
    }
    close(fd);
    data_size += len;
    buf.push_back(std::vector<int>());
    count.push_back(0);
  }

  template <class key_type, class value_type>
  typename MyKV<key_type, value_type>::MyBuf MyKV<key_type, value_type>::getbuf(size_t file_idx, size_t cs) {
    return MyBuf(this, file_idx, cs);
  }

  template <class key_type, class value_type>
  std::pair<const int*, size_t> MyKV<key_type, value_type>::getbuf(size_t idx) {
    return { buf[idx].data(), buf[idx].size() };
  }

  template <class key_type, class value_type>
  void MyKV<key_type, value_type>::getbuf(std::vector<int>& tmpbuf, size_t idx, size_t size) {
    if (size == 0) {
      // when size == 0, all the contents in file/idx are loaded to tmpbuf
      std::string fname;

      if (file_exist[idx] == 1) {
        // file exists
        int fd = open(fname.c_str(), O_RDONLY, (mode_t)0600);
        if (fd == -1) {
          perror("Error opening file");
          exit(EXIT_FAILURE);
        }
        off_t fsize = lseek(fd, (size_t)0, SEEK_END);
        lseek(fd, (size_t)0, SEEK_SET);
        size_t length = fsize / sizeof(int);
        //    assert(fsize == length * sizeof(int));
        size_t os = tmpbuf.size();
        tmpbuf.resize(tmpbuf.size() + length);
        ssize_t r = read(fd, tmpbuf.data() + os, fsize);
        if (r < 0) {
          perror("db error 1");
          exit(-1);
        }
        close(fd);
      }
      tmpbuf.insert(tmpbuf.end(), buf[idx].begin(), buf[idx].end());
      /* if (tmpbuf[0] != 0 && tmpbuf[0] != 1) {
        return;
      }*/
    }
    else {
      if (size > buf[idx].size()) {
        // when the required size is greater than the buffer, read the remaining from the file
        std::string fname = DBPath + "/" + std::to_string(idx) + ".dat";
        int fd = open(fname.c_str(), O_RDONLY, (mode_t)0600);
        if (fd == -1) {
          perror("Error opening file");
          exit(EXIT_FAILURE);
        }
        size_t os = tmpbuf.size();
        tmpbuf.resize(tmpbuf.size() + size);
        ssize_t r = read(fd, tmpbuf.data() + os, size * sizeof(int));
        if (r < 0) {
          perror("db error 1");
          exit(-1);
        }
        close(fd);
      }
      else {
        // when the required size is smaller than the buffer, the data can be directly read from the buffer
        tmpbuf.insert(tmpbuf.end(), buf[idx].begin(), buf[idx].begin() + size);
      }
    }
  }

  template <class key_type, class value_type>
  void MyKV<key_type, value_type>::merge(const key_type& k, const void* a, size_t len,
    bool store_value, double mni, const std::vector<std::vector<unsigned>>& orbits, const std::vector<unsigned>& perm, bool adaptive_sampling) {
    size_t file_id;
    bool first = false;
    //std::vector<int> rperm(len / sizeof(value_type) - 1);

    // check if the key is already in the KV_Store
    auto it = keys.find(k);
    if (it != keys.end()) {
      file_id = it->second;
    }
    else {
      file_id = nfiles++;
      keys[k] = file_id;
      file_exist.push_back(0);
      buf.emplace_back((int*)a, (int*)a + len / sizeof(value_type));
      if (mni >= 0) {

        std::vector<std::set<int>> vs;

        for (int i = 1; i < len / sizeof(value_type); i++) vs.push_back({ *((int*)a + i) });

        std::vector<std::set<int>> vos;
        for (auto& o : orbits) {
          std::set<int> vo;
          for (unsigned s : o) {
            vo.insert(vs[s].begin(), vs[s].end());
          }
          vos.push_back(vo);
        }

        std::vector<std::set<int>> vs_new(vs.size());


        for (int i = 0; i < vos.size(); i++) {
          for (unsigned ss : orbits[i]) {
            vs_new[perm[ss]] = vos[i];
          }
        }

        for (int si = 0; si < vs.size(); si++) {
          if (vs_new[perm[si]].size() == 0) {
            vs_new[perm[si]] = vs[si];
          }
        }

        if (adaptive_sampling)
          qp_set.push_back({ {*((int*)a), 1} });

        mni_met.push_back(false);
        distinct_vertices.push_back(vs_new);
      }
      count.push_back(1);

      first = true;
      data_size += len;
    }


    if (!first) {
      if (mni >= 0) {
        assert(ncols == len / sizeof(value_type));
        if (adaptive_sampling)
        {
          if (qp_set[file_id].find(*((int*)a)) == qp_set[file_id].end()) qp_set[file_id][*((int*)a)] = 0;
          qp_set[file_id][*((int*)a)] += 1;
        }

        if (!mni_met[file_id]) {

          std::vector<std::set<int>> vs;
          //if (perm.size() > 0) {
           // for (int i = 1; i < ncols; i++) vs.push_back({ *((int*)a + rperm[i - 1] + 1) });
         // }
          //else {
          for (int i = 1; i < ncols; i++) vs.push_back({ *((int*)a + i) });
          //}

          std::vector<std::set<int>> vos;
          for (auto& o : orbits) {
            std::set<int> vo;
            for (unsigned s : o) {
              vo.insert(vs[s].begin(), vs[s].end());
            }
            vos.push_back(vo);
          }

          std::vector<std::set<int>> vs_new(vs.size());

          for (int i = 0; i < vos.size(); i++) {
            for (unsigned ss : orbits[i]) {
              vs_new[perm[ss]] = vos[i];
            }
          }

          for (int si = 0; si < vs.size(); si++) {
            if (vs_new[perm[si]].size() == 0) {
              vs_new[perm[si]] = vs[si];
            }
          }


          auto& vss = distinct_vertices[file_id];
          bool tflag = true;
          for (int i = 0; i < ncols - 1; i++) {
            vss[i].insert(vs_new[i].begin(), vs_new[i].end());
            if (vss[i].size() < mni) {
              //mni_buf_size += sizeof(int);
              tflag = false;
            }
          }
          if (tflag) {
            mni_met[file_id] = true;
            vss.clear();
          }
        }
      }

      if (store_value) {
        for (int i = 0; i < len / sizeof(value_type); i++) {
          buf[file_id].push_back(*((int*)a + i));
        }
        data_size += len;
        //assert(buf[file_id].size() % 5 == 0);
        if (buf[file_id].size() >= BUF_SIZE * 1024l) {
          std::ofstream fd((DBPath + "/" + std::to_string(file_id) + ".dat").c_str(), std::ios::out | std::ios::binary | std::ios::ate);
          fd.write((char*)buf[file_id].data(), buf[file_id].size() * sizeof(value_type));
          fd.close();
          data_size -= buf[file_id].size() * sizeof(value_type);
          buf[file_id].clear();
          file_exist[file_id] = 1;
        }
      }

      count[file_id] += 1;
    }

    if (data_size > BUF_SIZE * 1024l) {
      for (int i = 0; i < nfiles; i++) {
        std::ofstream fd((DBPath + "/" + std::to_string(i) + ".dat").c_str(),
          std::ios::out | std::ios::binary | std::ios::ate);
        //assert(buf[i].size() % 5 == 0);
        fd.write((char*)buf[i].data(), buf[i].size() * sizeof(value_type));
        fd.close();
        buf[i].clear();
        file_exist[i] = 1;
        //std::vector<int>().swap(buf[i]);
      }
      data_size = 0;
    }
  }


  template <class key_type, class value_type>
  void MyKV<key_type, value_type>::merge_count(const key_type& k, size_t c) {
    auto it = keys.find(k);
    if (it != keys.end()) {
      size_t file_id = it->second;
      count[file_id] += c;
    }
    else {
      keys[k] = nfiles++;
      count.push_back(c);
    }
  }



  template <class key_type, class value_type>
  size_t MyKV<key_type, value_type>::size() {
    return keys.size();
  }

  template <class key_type, class value_type>
  size_t MyKV<key_type, value_type>::num_subgraphs() {
    size_t tot_num = 0;
    for (auto& k : keys) {
      tot_num += buf[k.second].size() / ncols;
    }
    return tot_num;
  }

  template <class key_type, class value_type>
  MyKV<key_type, value_type>::~MyKV() {
    clear();
  }

  template <class key_type, class value_type>
  void MyKV<key_type, value_type>::clear() {
    keys.clear();
    buf.clear();
    count.clear();
    mni_met.clear();
    distinct_vertices.clear();
    file_exist.clear();
    qp_set.clear();
    qp_path.clear();
    int s = system(("rm -rf " + DBPath).c_str());
  }

  template class MyKV<int>;
  template class MyKV<std::string>;
}  // namespace euler::db
