#pragma once
#ifndef FILE_HANDLER
#define FILE_HANDLER
#include <filesystem>
#include <stdexcept>
#include <string>

namespace cactus {

class StructuresDir {
  private:
    std::string path;

    void validateDir(const std::string& dirPath) {
        if (!std::filesystem::exists(dirPath) || !std::filesystem::is_directory(dirPath))
            throw std::runtime_error("Invalid directory: " + dirPath);
    }

  public:
    StructuresDir(const std::string& dirPath) : path(dirPath) {
        validateDir(dirPath);
    }

    std::filesystem::directory_iterator dirFiles() {
        return std::filesystem::directory_iterator(path);
    }
  };
} // namespace cactus
#endif
