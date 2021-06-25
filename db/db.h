#pragma once

#include <assert.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <atomic>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace euler::db {

  const size_t BUF_SIZE = 1024l * 1024l * 128l;

  static std::vector<std::vector<unsigned>> dummy1 = {};
  static std::vector<unsigned> dummy2 = {};

  template <class key_type, class value_type = int>
  struct MyKV {
    static int db_id;
    std::string DBPath;
    std::atomic<size_t> nfiles;
    std::map<key_type, size_t> keys;
    std::vector<std::vector<value_type>> buf;
    std::vector<double> count;
    std::vector<bool> mni_met;
    int ncols;
    size_t data_size = 0; //, mni_buf_size = 0;
    std::vector<std::vector<std::set<value_type>>> distinct_vertices;
    std::vector<int> file_exist;

    std::vector<std::set<int>> qp_set;
    std::vector<std::vector<std::vector<int>>> qp_path;

    void print() {
      for (auto& [key, value] : keys) {
        std::cout << "key: " << key << std::endl;
        for (int i = 0; i < buf[value].size(); i++) {
          std::cout << buf[value][i] << " ";
          if ((i + 1) % ncols == 0) std::cout << std::endl;
        }
      }
    }

    void get_pattern_path(const std::vector<std::vector<int>>& qp_count) {

      for (auto& s : qp_set) {
        std::vector<std::vector<int>> res(s.size());
        int t = 0;
        for (auto ss: s) {
          int cur = ss;
          for (int i = qp_count.size() - 1; i > 0; i--) {
            int left = qp_count[i][2 * cur];
            int right = qp_count[i][2 * cur + 1];
            res[t].push_back(right);
            cur = left;
          }
          int left = qp_count[0][2 * cur];
          int right = qp_count[0][2 * cur + 1];
          res[t].push_back(right);
          res[t].push_back(left);
          t++;
        }
        qp_path.push_back(res);
      }
    }

    void combine(MyKV& other, bool mni, bool store, bool adaptive_sampling) {
      for (auto& [key, value] : other.keys) {
        if (keys.find(key) == keys.end()) {
          buf.push_back(std::vector<value_type>());
          for (auto& b : other.buf[value]) {
            buf[buf.size() - 1].push_back(b);
          }
          count.push_back(other.count[value]);
          keys[key] = nfiles;
          nfiles++;

          if (mni) {
            mni_met.push_back(other.mni_met[value]);
            distinct_vertices.push_back(other.distinct_vertices[value]);

            if (adaptive_sampling) {
              qp_path.push_back(other.qp_path[value]);
            }
          }
        }
        else {
          size_t fid = keys[key];
          if (store) {
            for (auto& b : other.buf[value]) {
              buf[fid].push_back(b);
            }
          }
          count[fid] += other.count[value];
          if (mni) {
            mni_met[fid] = (mni_met[fid] || other.mni_met[value]);
            if (!mni_met[fid]) {
              for (int i = 0; i < distinct_vertices[fid].size(); i++) {
                distinct_vertices[fid][i].insert(other.distinct_vertices[value][i].begin(), other.distinct_vertices[value][i].end());
              }
            }
            if (adaptive_sampling) {
              qp_path[fid].insert(qp_path[fid].end(), other.qp_path[value].begin(), other.qp_path[value].end());
            }
            //qp_set[fid].insert(other.qp_set[value].begin(), other.qp_set[value].end());
          }
        }
      }
      other.keys.clear();
      other.buf.clear();
      other.count.clear();
      other.distinct_vertices.clear();
    }

    // this is constant view buf and does not make any changes to db
    class MyBuf {
      MyKV* db;
      size_t file_idx;
      size_t total_length = 0;
      size_t chunk_size;

    public:
      MyBuf(MyKV* _db, size_t fidx, size_t cs) : db(_db), file_idx(fidx), chunk_size(cs) {
        size();
      }

      struct iterator {
        size_t MAX_BUF_SIZE;
        size_t local_offset, global_offset;
        int* buffer;
        size_t buffer_size;
        MyBuf& mybuf;
        bool in_file = false;  // indicate whether the current data is read from file
        bool has_next = true;

        iterator(MyBuf& b) : mybuf(b), buffer(nullptr), buffer_size(0), local_offset(0), global_offset(0), MAX_BUF_SIZE(0) {}

        ~iterator() {
          // No delete is needed, because the last buffer is a pointer to a vector
          // delete[] buffer;
        }

        iterator& operator++() {
          if (local_offset < buffer_size) {
            local_offset += mybuf.chunk_size;
            global_offset += mybuf.chunk_size;
          }
          if (local_offset == buffer_size && global_offset != mybuf.total_length) {
            local_offset = 0;
            buffer_size = 0;
            if (in_file) {
              std::string fname = mybuf.db->DBPath + "/" + std::to_string(mybuf.file_idx) + ".dat";
              // file exists
              int fd = open(fname.c_str(), O_RDONLY, (mode_t)0600);
              if (fd == -1) {
                perror("Error opening file");
                exit(EXIT_FAILURE);
              }
              size_t read_size = MAX_BUF_SIZE > (mybuf.total_length - global_offset) ? (mybuf.total_length - global_offset) : MAX_BUF_SIZE;
              lseek(fd, global_offset * sizeof(value_type), SEEK_SET);
              buffer_size = read(fd, buffer, read_size * sizeof(value_type));
              buffer_size /= sizeof(value_type);
              //   if (buffer[0] != 0 && buffer[0] != 1) {
              //    std::cout << "wwww" << std::endl;
              // }
              if (buffer_size == 0) {
                if (in_file) {
                  delete buffer;
                  in_file = false;
                }
                buffer = mybuf.db->buf[mybuf.file_idx].data();
                MAX_BUF_SIZE = buffer_size = mybuf.db->buf[mybuf.file_idx].size();
              }
              close(fd);
            }
          }
          return *this;
        }

        // load next buffer
        void next() {
          global_offset += buffer_size;
          //std::cout << "global offset: " << global_offset << std::endl;
          local_offset = 0;
          buffer_size = 0;
          if (in_file) {
            std::string fname = mybuf.db->DBPath + "/" + std::to_string(mybuf.file_idx) + ".dat";
            // file exists
            // TODO: some bug here
            int fd = open(fname.c_str(), O_RDONLY, (mode_t)0600);
            if (fd == -1) {
              perror("Error opening file");
              exit(EXIT_FAILURE);
            }
            size_t read_size = MAX_BUF_SIZE > (mybuf.total_length - global_offset) ? (mybuf.total_length - global_offset) : MAX_BUF_SIZE;
            lseek(fd, global_offset * sizeof(value_type), SEEK_SET);
            buffer_size = read(fd, buffer, read_size * sizeof(value_type));
            buffer_size /= sizeof(value_type);
            //   if (buffer[0] != 0 && buffer[0] != 1) {
            //    std::cout << "wwww" << std::endl;
            // }
            close(fd);
            if (buffer_size == 0) {
              if (in_file) {
                delete buffer;
                in_file = false;
              }
              buffer = mybuf.db->buf[mybuf.file_idx].data();
              MAX_BUF_SIZE = buffer_size = mybuf.db->buf[mybuf.file_idx].size();
              has_next = false;
            }
          }
        }

        bool operator!=(const iterator& it) const {
          return global_offset != it.global_offset;
        }
        bool operator==(const iterator& it) const {
          return global_offset == it.global_offset;
        }

        // return the starting pointer
        value_type* operator*() {
          return buffer + local_offset;
        }
      };

    public:
      iterator begin() {
        iterator it(*this);
        std::string fname;

        if (db->file_exist[file_idx]) {
          std::string fname = db->DBPath + "/" + std::to_string(file_idx) + ".dat";
          // file exists
          it.in_file = true;
          it.MAX_BUF_SIZE = BUF_SIZE / chunk_size * chunk_size;
          it.buffer = new value_type[it.MAX_BUF_SIZE];
          it.buffer_size = 0;

          int fd = open(fname.c_str(), O_RDONLY, (mode_t)0600);
          if (fd == -1) {
            perror("Error opening file");
            exit(EXIT_FAILURE);
          }
          it.local_offset = it.global_offset = 0;
          size_t read_size = it.MAX_BUF_SIZE > total_length ? total_length : it.MAX_BUF_SIZE;
          size_t buffer_size = read(fd, it.buffer, read_size * sizeof(value_type));
          if (buffer_size > 0) {
            it.buffer_size = buffer_size / sizeof(value_type);
          }
          close(fd);
        }
        else {
          it.in_file = false;
          it.buffer = db->buf[file_idx].data();
          it.MAX_BUF_SIZE = it.buffer_size = db->buf[file_idx].size();
          it.has_next = false;
        }

        return it;
      }

      iterator end() {
        iterator it(*this);
        it.global_offset = total_length;
        return it;
      }

      size_t size() {
        if (total_length == 0) {
          total_length = db->buf[file_idx].size();

          if (db->file_exist[file_idx]) {
            // file exists
            std::string fname = db->DBPath + "/" + std::to_string(file_idx) + ".dat";
            int fd = open(fname.c_str(), O_RDONLY, (mode_t)0600);
            if (fd == -1) {
              perror("Error opening file");
              exit(EXIT_FAILURE);
            }
            off_t fsize = lseek(fd, (size_t)0, SEEK_END);
            //std::cout << "fsize: " << fsize << std::endl;
            lseek(fd, (size_t)0, SEEK_SET);
            total_length += fsize / sizeof(value_type);
            close(fd);
          }
          assert(total_length % chunk_size == 0);
        }
        return total_length;
      }
    };

  public:
    MyKV(int nc);

    void put(const void* a, size_t len);

    void merge(const key_type& k, const void* a, size_t len,
      bool store_value = true, size_t mni = 0, const std::vector<std::vector<unsigned>>& orbits = dummy1, const std::vector<unsigned>& perm = dummy2);

    MyBuf getbuf(size_t file_idx, size_t cs);

    void getbuf(std::vector<int>& tmpbuf, size_t idx, size_t size = 0);

    std::pair<const int*, size_t> getbuf(size_t idx);

    //void replace(const key_type &k, const void *a, size_t len);

    size_t size();

    ~MyKV();
  };
}  // namespace euler::db
