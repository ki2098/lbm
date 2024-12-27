#include <cstdio>

int main() {
    int sum = 100;
    #pragma acc parallel loop reduction(+:sum)
    for (int i = 0; i < 5000; i ++) {
        sum -= 1;
    }
    printf("%d\n", sum);
}