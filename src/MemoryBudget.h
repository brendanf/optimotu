#ifndef OPTIMOTU_MEMORYBUDGET_H_INCLUDED
#define OPTIMOTU_MEMORYBUDGET_H_INCLUDED

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <string>

class MemoryBudgetExceeded : public std::runtime_error {
public:
  const std::size_t budget_bytes;
  const std::size_t current_bytes;
  const std::size_t requested_bytes;
  const std::string context;

  MemoryBudgetExceeded(
    std::size_t budget_bytes,
    std::size_t current_bytes,
    std::size_t requested_bytes,
    const std::string &context
  ) :
    std::runtime_error("clustering memory budget exceeded"),
    budget_bytes(budget_bytes),
    current_bytes(current_bytes),
    requested_bytes(requested_bytes),
    context(context) {}
};

class MemoryBudgetTracker {
private:
  const std::size_t budget_bytes;
  std::size_t current_bytes = 0;
  std::size_t peak_bytes = 0;
  mutable std::mutex mutex;

public:
  explicit MemoryBudgetTracker(std::size_t budget_bytes) :
    budget_bytes(budget_bytes) {}

  bool enabled() const {
    return budget_bytes > 0;
  }

  std::size_t budget() const {
    return budget_bytes;
  }

  std::size_t current() const {
    std::lock_guard<std::mutex> lock(mutex);
    return current_bytes;
  }

  std::size_t peak() const {
    std::lock_guard<std::mutex> lock(mutex);
    return peak_bytes;
  }

  void acquire(std::size_t bytes, const std::string &context) {
    if (!enabled() || bytes == 0) {
      return;
    }
    std::lock_guard<std::mutex> lock(mutex);
    if (current_bytes > budget_bytes || bytes > (budget_bytes - current_bytes)) {
      throw MemoryBudgetExceeded(budget_bytes, current_bytes, bytes, context);
    }
    current_bytes += bytes;
    if (current_bytes > peak_bytes) {
      peak_bytes = current_bytes;
    }
  }

  void release(std::size_t bytes) {
    if (!enabled() || bytes == 0) {
      return;
    }
    std::lock_guard<std::mutex> lock(mutex);
    if (bytes >= current_bytes) {
      current_bytes = 0;
    } else {
      current_bytes -= bytes;
    }
  }
};

#endif // OPTIMOTU_MEMORYBUDGET_H_INCLUDED
