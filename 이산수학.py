import copy

# ===================================================================
# 1. 행렬 입력 기능
# ===================================================================
def get_matrix_input():
    """사용자로부터 n x n 행렬을 입력받아 반환합니다."""
    while True:
        try:
            n = int(input("행렬의 크기 n을 입력하세요 (예: 3): "))
            if n <= 0:
                print("n은 0보다 큰 정수여야 합니다.")
                continue
            break
        except ValueError:
            print("숫자를 입력해주세요.")  # 잘못된 입력시 다시 반복

    print(f"{n}x{n} 행렬의 원소를 한 행씩 입력하세요 (각 원소는 공백으로 구분).")
    matrix = []
    for i in range(n):
        while True:
            try:
                row_input = input(f"{i + 1}번째 행 입력: ")
                row = [float(x) for x in row_input.split()]  # n개의 원소를 받아 행렬저장
                if len(row) != n:
                    print(f"정확히 {n}개의 원소를 입력해야 합니다.")
                    continue
                matrix.append(row)
                break
            except ValueError:
                print("원소는 숫자여야 합니다. 다시 입력해주세요.")  # 잘못된 입력시 다시 반복
    return matrix

# ===================================================================
# 2. 행렬식을 이용한 역행렬 계산 기능
# ===================================================================
def get_determinant(matrix):
    """재귀 함수를 이용해 행렬식을 계산합니다 (여인수 전개 사용)."""
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    
    # 여인수 전개
    determinant = 0
    for j in range(n):
        sign = (-1) ** j
        # 1행을 기준으로 소행렬 구하기, j열 제외하고 구함
        minor = [row[:j] + row[j+1:] for row in matrix[1:]]
        determinant += sign * matrix[0][j] * get_determinant(minor)  # 소행렬식 계산
    return determinant

def get_adjugate_matrix(matrix):
    """수반 행렬(Adjugate Matrix)을 계산합니다."""
    n = len(matrix)
    if n == 1:
        # 0x0 행렬의 행렬식은 1이므로 크기가 1일때의 수반행렬은 항상 1
        return [[1]]

    adjugate = []
    for i in range(n):
        cofactor_row = []
        for j in range(n):
            sign = (-1) ** (i + j)
            # 소행렬식 계산
            minor = [row[:j] + row[j+1:] for row in (matrix[:i] + matrix[i+1:])]
            cofactor = sign * get_determinant(minor)  # 여인수 행렬의 원소 계산
            cofactor_row.append(cofactor)  # 여인수 행렬에 추가
        adjugate.append(cofactor_row)

    # 여인수 행렬을 전치해 수반행렬로 변환
    adjugate_transposed = [[adjugate[j][i] for j in range(n)] for i in range(n)]
    return adjugate_transposed

def inverse_by_determinant(matrix):
    """행렬식과 수반 행렬을 이용해 역행렬을 구합니다."""
    # 행렬식을 통해 주어진 행렬이 역행렬이 존재하는지 확인후 0이 아니라면 역행렬 리턴
    determinant = get_determinant(matrix)
    if determinant == 0:
        return None

    adjugate = get_adjugate_matrix(matrix)
    n = len(matrix)
    # 역행렬에 (수반행렬/행렬식)의 원소 채워넣기
    inverse = [[adjugate[i][j] / determinant for j in range(n)] for i in range(n)]
    return inverse

# ===================================================================
# 3. 가우스-조던 소거법을 이용한 역행렬 계산 기능
# ===================================================================
def inverse_by_gauss_jordan(matrix):
    """가우스-조던 소거법을 이용해 역행렬을 구합니다."""
    n = len(matrix)
    # 기존 행렬 보존하기위해 딥카피 수행
    mat = copy.deepcopy(matrix)
    # i=j(주대각 성분)일 때 부울함수를 float()를 통해 1로 만들어 단위행렬 생성
    identity = [[float(i == j) for j in range(n)] for i in range(n)]

    for i in range(n):
        # 피벗(주대각 성분으로 설정)이 0일 경우, 아래 행과 교환하여 0이 아니게 만듦
        if mat[i][i] == 0:
            for j in range(i + 1, n):
                if mat[j][i] != 0:
                    mat[i], mat[j] = mat[j], mat[i]
                    identity[i], identity[j] = identity[j], identity[i] # 단위행렬에도 같은 연산수행
                    break
            else: # 교환할 행이 없으면 역행렬이 존재하지 않음(열이 000이라는 소리이므로)
                return None

        # 피벗을 1로 만들기
        pivot = mat[i][i]
        for j in range(i, n):
            mat[i][j] /= pivot
        for j in range(n):
            identity[i][j] /= pivot

        # 현재 열의 다른 원소들을 0으로 만들기
        for j in range(n):
            if i != j:
                factor = mat[j][i]
                for k in range(i, n):
                    mat[j][k] -= factor * mat[i][k]#기본행 연산으로 factor를 0으로 만듦
                for k in range(n):
                    identity[j][k] -= factor * identity[i][k] # 단위행렬에도 같은연산수행

    return identity # 가우스 조르단 소거법으로 구해진 역행렬 리턴

# ===================================================================
# 4.  출력, 비교 기능
# ===================================================================
def print_matrix(matrix, name="행렬"):
    """행렬을 지정된 이름과 함께 소수점 4자리까지 출력합니다."""
    if matrix is None:
        print(f"\n{name}은(는) 존재하지 않습니다 (행렬식이 0).")
        return
    print(f"\n--- {name} ---")
    for row in matrix:
        print(" ".join(f"{elem:10.4f}" for elem in row))

def compare_matrices(mat1, mat2, tolerance=1e-9):
    """두 행렬이 오차 범위(tolerance) 내에서 동일한지 비교합니다.{(1/3)같은 경우에 값이 다르다고 나올확률 있으므로 오차범위(1e-9) 설정}"""
    if mat1 is None or mat2 is None:
        return mat1 is None and mat2 is None
    if len(mat1) != len(mat2) or len(mat1[0]) != len(mat2[0]):
        return False
    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            if abs(mat1[i][j] - mat2[i][j]) > tolerance:
                return False
    return True

# ===================================================================
# 5. 추가 기능 (검산)
# ===================================================================
def multiply_matrices(mat1, mat2):
    """두 행렬을 곱한 결과를 반환합니다."""
    n = len(mat1)
    result = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                result[i][j] += mat1[i][k] * mat2[k][j]
    return result

def verify_inverse(original_matrix, inverse_matrix):
    """원본 행렬과 역행렬을 곱하여 단위 행렬이 나오는지 검산합니다."""
    if inverse_matrix is None:
        return None # 역행렬이 없으므로 검산 불가

    n = len(original_matrix)
    
    # 원본 행렬과 역행렬을 곱함
    product = multiply_matrices(original_matrix, inverse_matrix)

    # 단위 행렬 생성
    identity = [[float(i == j) for j in range(n)] for i in range(n)]

    # 곱셈 결과가 단위 행렬과 같은지 비교
    is_correct = compare_matrices(product, identity)

    # 검산 과정을 보여주기 위해 곱셈 결과도 함께 반환
    return is_correct, product

# ===================================================================
# 6. 메인 실행 함수
# ===================================================================
def main():
    """프로그램의 주 실행 로직을 담당합니다."""
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

    # 검산 결과 출력
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


