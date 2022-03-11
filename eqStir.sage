### Check equivariant log concavity of stirling numbers
### Authors: Siddarth Kannan, Shiyue Li

import warnings

def frob_char_ass(n):
    p = SymmetricFunctions(QQ).powersum()
    s = SymmetricFunctions(QQ).schur()
    div = divisors(n)
    codiv = []
    frob = sum(euler_phi(d)/n * (p[d])^(n//d) for d in div)
    return s(frob)


def frob_first(n,m):
    h = SymmetricFunctions(QQ).homogeneous()
    s = SymmetricFunctions(QQ).schur()
    insert = sum(frob_char_ass(i) for i in range(1, n - m + 2))
    total = h[m].plethysm(insert)
    ans = total.restrict_degree(n)
    return s(ans)

def frob_second(n,m):
#this calculates the character of the S_n rep spanned by m-partitions of {1,...,n}, as a symmetric function
#the coefficient of s_lambda is the multiplicity of the corresponding irrep
    Sym = SymmetricFunctions(QQ)
    s = SymmetricFunctions(QQ).schur()
    h = Sym.homogeneous()
    src = h[m]
    insert = sum(h[i] for i in range(1, n - m + 2))
    total = src.plethysm(insert)
    ans= total.restrict_degree(n)
    return s(ans)

def is_subrep(f, g, n):
#f and g are symmetric functions, homogeneous of degree n, representing the characters of V and W, respectively
#checks if V is an S_n-subrep of W
    Sym = SymmetricFunctions(QQ)
    s = SymmetricFunctions(QQ).schur()
    potential_subrep = s(f)
    potential_contain = s(g)
    partition_list = Partitions(n)
    for i in range(0, len(partition_list)):
        if potential_subrep.scalar(s[partition_list[i]]) > potential_contain.scalar(s[partition_list[i]]):
            print(i, "-th term, subrep coeff: ", potential_subrep.scalar(s[partition_list[i]]))
            print(i, "-th term, bigrep coeff: ", potential_contain.scalar(s[partition_list[i]]))
            print(partition_list[i])
            return False
    return True


def dim(f,n):
    #computes the dimension of a S_n representation using the Frobenius
        p = SymmetricFunctions(QQ).powersum()
        part = []
        for i in range(1, n+1):
            part.append(1)
        return f.scalar(p(part))

def eq_first(n):
#checks if the graded v.s. spanned by partitions of {1,...,n} is equiv log concave
#grading is given by number of blocks in partition
    for i in range(2, n):
        left = frob_first(n, i-1)
        right = frob_first(n, i+1)
        middle = frob_first(n, i)

        #check dimensions
        dim_left = dim(left, n)
        dim_right = dim(right, n)
        dim_middle = dim(middle, n)
        print(i,":dim left is", dim_left, ", dim right is", dim_right, ", dim middle is", dim_middle)
        if (dim_left != stirling_number1(n, i-1) or
            dim_right != stirling_number1(n, i+1) or
            dim_middle != stirling_number1(n, i)):
            warnings.warn("Dimensions are wrong for i = "+str(i))

        test = left.itensor(right)
        check = middle.itensor(middle)

        if not is_subrep(test, check, n):
            print("Not equiv log concave at i = "+str(i))
    return True


def eq_second(n):
#checks if the graded v.s. spanned by partitions of {1,...,n} is equiv log concave
#grading is given by number of blocks in partition
      for i in range(2, n):
          left = frob_second(n, i-1)
          right = frob_second(n, i+1)
          middle = frob_second(n, i)

          #check dimensions
          dim_left = dim(left, n)
          dim_right = dim(right, n)
          dim_middle = dim(middle, n)
          print(i,":dim left is", dim_left, ", dim right is", dim_right, ", dim middle is", dim_middle)
          if (dim_left != stirling_number2(n, i-1) or
              dim_right != stirling_number2(n, i+1) or
              dim_middle != stirling_number2(n, i)):
              warnings.warn("Dimensions are wrong for i = "+str(i))

          test = left.itensor(right)
          check = middle.itensor(middle)

          if not is_subrep(test, check, n):
              print("Not equiv log concave at i = "+str(i))

      return True
