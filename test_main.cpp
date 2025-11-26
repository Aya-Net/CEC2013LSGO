//
// Created by qq328 on 2025/11/26.
//
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <atomic>
#include <iostream>

static bool preallocateFilePosix(const std::string &path, uint64_t size_bytes, int &outFd) {
  int fd = open(path.c_str(), O_CREAT | O_RDWR, 0644);
  if (fd < 0) return false;
  int rc = posix_fallocate(fd, 0, (off_t)size_bytes);
  if (rc != 0) {
    // 回退：尝试 ftruncate（注意：fallocate 可能在不同系统上失败）
    if (ftruncate(fd, (off_t)size_bytes) != 0) {
      close(fd);
      return false;
    }
  }
  outFd = fd;
  return true;
}
int main() {
  std::string pop_file_path = "population.bin";

  const uint64_t slot_bytes = (uint64_t)sizeof(double) * (uint64_t)1000;
  const uint64_t PREALLOC_BYTES = 600ULL * 1024ULL * 1024ULL * 1024ULL; // 600GB
  uint64_t capacity = PREALLOC_BYTES / slot_bytes;
  if (capacity == 0) capacity = 1;

  int pop_fd = -1;
  if (!preallocateFilePosix(pop_file_path, capacity * slot_bytes, pop_fd)) {
    std::cerr << "Failed to create/preallocate `population.bin`\n";
    return 0;
  }
}