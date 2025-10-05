import copy


# 1. 행렬 입력 기능
def get_matrix_input():
    while True:
        try:
            n = int(input("행렬의 크기 n을 입력하세요 (예: 3): "))
            if n <= 0:
                print("n은 0보다 큰 정수여야 합니다.")
                continue
            break
        except ValueError:
            print("숫자를 입력해주세요.")#잘못된 입력시 다시 반복

    print(f"{n}x{n} 행렬의 원소를 한 행씩 입력하세요 (각 원소는 공백으로 구분).")
    matrix = []
    for i in range(n):
        while True:
            try:
                row_input = input(f"{i + 1}번째 행 입력: ")
                row = [float(x) for x in row_input.split()]#n개의 원소를 받아 행렬저장
                if len(row) != n:
                    print(f"정확히 {n}개의 원소를 입력해야 합니다.")
                    continue
                matrix.append(row)
                break
            except ValueError:
                print("원소는 숫자여야 합니다. 다시 입력해주세요.")#잘못된 입력시 다시 반복
    return matrix

# 2. 행렬식을 이용한 역행렬 계산 기능

def get_determinant(matrix):
    n = len(matrix)
    if n == 1:
        return matrix[0][0]

    determinant = 0
    for j in range(n):
        sign = (-1) ** j
        minor = [row[:j] + row[j+1:] for row in matrix[1:]]# 소행렬 구하기
        determinant += sign * matrix[0][j] * get_determinant(minor)
    return determinant

def get_adjugate_matrix(matrix):
    n = len(matrix)
    if n == 1:
        return [[1]]

    adjugate = []
    for i in range(n):
        cofactor_row = []
        for j in range(n):
            sign = (-1) ** (i + j)
            minor = [row[:j] + row[j+1:] for row in (matrix[:i] + matrix[i+1:])]
            cofactor = sign * get_determinant(minor)
            cofactor_row.append(cofactor)
        adjugate.append(cofactor_row)

    adjugate_transposed = [[adjugate[j][i] for j in range(n)] for i in range(n)]
    return adjugate_transposed

def inverse_by_determinant(matrix):
    determinant = get_determinant(matrix)
    if determinant == 0:
        return None

    adjugate = get_adjugate_matrix(matrix)
    n = len(matrix)
    inverse = [[adjugate[i][j] / determinant for j in range(n)] for i in range(n)]
    return inverse

 
# 3. 가우스-조던 소거법을 이용한 역행렬 계산 기능
 
def inverse_by_gauss_jordan(matrix):
    n = len(matrix)
    mat = copy.deepcopy(matrix)
    identity = [[float(i == j) for j in range(n)] for i in range(n)]

    for i in range(n):
        if mat[i][i] == 0:
            for j in range(i + 1, n):
                if mat[j][i] != 0:
                    mat[i], mat[j] = mat[j], mat[i]
                    identity[i], identity[j] = identity[j], identity[i]
                    break
            else:
                return None

        pivot = mat[i][i]
        for j in range(i, n):
            mat[i][j] /= pivot
        for j in range(n):
            identity[i][j] /= pivot

        for j in range(n):
            if i != j:
                factor = mat[j][i]
                for k in range(i, n):
                    mat[j][k] -= factor * mat[i][k]
                for k in range(n):
                    identity[j][k] -= factor * identity[i][k]

    return identity

 
# 4. 결과 출력 및 비교 기능
 
def print_matrix(matrix, name="행렬"):
    if matrix is None:
        print(f"\n{name}은(는) 존재하지 않습니다 (행렬식이 0).")
        return
    print(f"\n--- {name} ---")
    for row in matrix:
        print(" ".join(f"{elem:10.4f}" for elem in row))

def compare_matrices(mat1, mat2, tolerance=1e-9):
    if mat1 is None or mat2 is None:
        return mat1 is None and mat2 is None
    if len(mat1) != len(mat2) or len(mat1[0]) != len(mat2[0]):
        return False
    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            if abs(mat1[i][j] - mat2[i][j]) > tolerance:
                return False
    return True

 
# ✨ 추가 기능: 역행렬 검산
 
def multiply_matrices(mat1, mat2):
    n = len(mat1)
    result = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                result[i][j] += mat1[i][k] * mat2[k][j]
    return result

def verify_inverse(original_matrix, inverse_matrix):
    if inverse_matrix is None:
        return None # 검산 불가

    n = len(original_matrix)
    # 원본 행렬과 역행렬을 곱함
    product = multiply_matrices(original_matrix, inverse_matrix)

    # 단위 행렬 생성
    identity = [[float(i == j) for j in range(n)] for i in range(n)]

    # 곱셈 결과가 단위 행렬과 같은지 비교
    is_correct = compare_matrices(product, identity)

    # 검산 과정을 보여주기 위해 곱셈 결과도 함께 반환
    return is_correct, product

 
# 메인 실행 함수
 
def main():
    input_matrix = get_matrix_input()
    print_matrix(input_matrix, "입력된 행렬")

    inv_det = inverse_by_determinant(input_matrix)
    print_matrix(inv_det, "방법 1: 행렬식을 이용한 역행렬")

    inv_gj = inverse_by_gauss_jordan(input_matrix)
    print_matrix(inv_gj, "방법 2: 가우스-조던 소거법을 이용한 역행렬")

    print("\n--- 결과 비교 ---")
    if compare_matrices(inv_det, inv_gj):
        print("두 방법으로 계산된 역행렬의 결과가 동일합니다.")
    else:
        print("두 방법으로 계산된 역행렬의 결과가 다릅니다.")

    # ✨ 추가된 검산 결과 출력
    print("\n--- 역행렬 검산 (A * A^-1 = I) ---")
    verification_result = verify_inverse(input_matrix, inv_gj)
    if verification_result:
        is_correct, product_matrix = verification_result
        print_matrix(product_matrix, "원본 행렬 * 역행렬 결과")
        if is_correct:
            print("✅ 검산 결과: 정확한 역행렬입니다.")
        else:
            print("❌ 검산 결과: 역행렬 계산에 오류가 있습니다.")
    else:
        print("역행렬이 존재하지 않아 검산을 수행할 수 없습니다.")


if __name__ == "__main__":
    main()