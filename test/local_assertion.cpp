#include <local_assertion.hpp>

int main() {
    int x = 3;
    LOCAL_ASSERT_IN_RANGE(x, 3, 3);
    LOCAL_ASSERT_IN_RANGE(x, 0, 1);
}