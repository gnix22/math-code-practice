def collatz_function(n: int):
    print(f"Collatz Sequence for n={n}:")
    print(n)
    while n != 1:
        if n % 2 == 0:
            iter = n/2
        else:
            iter = 3*n+1
        n = iter
        print(iter)
