#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <vector>
#include <numeric>
using namespace std;

class Matrix {
	private:
		vector<double> data;
		int rows, cols;
	public:
		// Конструкторы
		Matrix() : rows(0), cols(0) {}
		Matrix(int r, int c, double val = 0.0) : rows(r), cols(c), data(r * c, val) {}
		// Доступ к элементам
		double& operator()(int i, int j) {
			return data[i * cols + j];
		}
		const double& operator()(int i, int j) const {
			return data[i * cols + j];
		}
		// Возвращение размеров
		int getRows() const {
			return rows;
		}
		int getCols() const {
			return cols;
		}
		// Сложение двух матриц
		Matrix operator+(const Matrix& other) const {
			if (rows != other.rows || cols != other.cols) // Проверка на одинаковую размерность матриц
				throw invalid_argument("Matrix dimensions do not match for addition");
			Matrix result(rows, cols);
			for (int i = 0; i < rows * cols; ++i)
				result.data[i] = data[i] + other.data[i];
			return result;
		}
		// Вычитание матриц
		Matrix operator-(const Matrix& other) const {
			if (rows != other.rows || cols != other.cols)// Проверка на одинаковую размерность матриц
				throw invalid_argument("Matrix dimensions do not match for subtraction");
			Matrix result(rows, cols);
			for (int i = 0; i < rows * cols; ++i)
				result.data[i] = data[i] - other.data[i];
			return result;
		}
		// Умножение матриц
		Matrix operator*(const Matrix& other) const {
			if (cols != other.rows)// Проверка на корректную размерность матриц
				throw invalid_argument("Matrix dimensions do not match for multiplication");
			Matrix result(rows, other.cols);
			for (int i = 0; i < rows; ++i)
				for (int k = 0; k < cols; ++k)
					for (int j = 0; j < other.cols; ++j)
						result(i, j) += (*this)(i, k) * other(k, j);
			return result;
		}
		// Умножение на скаляр
		Matrix operator*(double scalar) const {
			Matrix result(rows, cols);
			for (int i = 0; i < rows * cols; ++i)
				result.data[i] = data[i] * scalar;
			return result;
		}
		friend Matrix operator*(double scalar, const Matrix& m) {
			return m * scalar;
		}
		// Транспонирование
		Matrix transpose() const {
			Matrix result(cols, rows);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					result(j, i) = (*this)(i, j);
			return result;
		}
		// Чтение матрицы из файла
		void readFromFile(const string& filename) {
			ifstream file(filename);
			if (!file) throw runtime_error("Error opening file: " + filename);
			// Инициализация счетчиков строк и столбцов
			rows = 0;
			cols = 0;
			string line;
			vector<vector<double>> tempData;
			// Читаем файл построчно
			while (getline(file, line)) {
				if (line.empty()) continue;
				istringstream iss(line); // Используем строковый поток
				vector<double> row;
				double val;
				while (iss >> val) row.push_back(val);// Читаем все числа из строки
				if (rows == 0) cols = row.size(); // Если первая строка, то запоминаем размер
				else if ((int)row.size() != cols) // Проверка несогласованности количества столбцов
					throw runtime_error("Inconsistent number of columns in matrix");
				tempData.push_back(row);
				rows++;
			}
			// Проверка на пустоту
			if (rows == 0 || cols == 0)
				throw runtime_error("Empty matrix in file");
			// Изменяем размер хранилища данных и заполняем его
			data.resize(rows * cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					(*this)(i, j) = tempData[i][j];
		}
		// Вывод матрицы в поток
		friend ostream& operator<<(ostream& os, const Matrix& m) {
			int max_width = 0;
			for (int i = 0; i < m.rows; ++i)
				for (int j = 0; j < m.cols; ++j) {
					double val = abs(m(i, j)) < 1e-10 ? 0 : m(i, j); // Обнуляем малые числа
					ostringstream oss; // Строковый поток для определение ширины
					if (val == floor(val)) oss << (int)val; // Целые как целые выводим
					else oss << val;
					max_width = max(max_width, (int)oss.str().length()); // Обновляем максимальную ширину
				}
			// Выводим также матрицу с выравниванием
			for (int i = 0; i < m.rows; ++i) {
				for (int j = 0; j < m.cols; ++j) {
					double val = abs(m(i, j)) < 1e-10 ? 0 : m(i, j);
					if (val == floor(val)) os << setw(max_width + 1) << (int)val;
					else os << setw(max_width + 1) << val;
				}
				os << endl;
			}
			return os;
		}
		// Вычисление определителя
		double determinant() const {
			// Проверка на квадратность
			if (rows != cols) throw invalid_argument("Matrix must be square");
			Matrix LU(*this);
			double det = 1;
			int n = rows;
			for (int i = 0; i < n; ++i) {
				int maxRow = i;
				// Определяем строку с максимальным по модуля элементов в текущем столбце
				for (int k = i + 1; k < n; ++k)
					if (abs(LU(k, i)) > abs(LU(maxRow, i))) maxRow = k;
				if (maxRow != i) {// Перестановка строк
					swap(LU(i, i), LU(maxRow, i));
					det *= -1;
				}
				if (abs(LU(i, i)) < 1e-12) return 0;// Проверка на вырожденность
				det *= LU(i, i); // Умножение определителя на диогональный элемент
				for (int k = i + 1; k < n; ++k) {
					LU(k,i) /= LU(i,i); // Сохраняем множитель
					for (int j = i + 1; j < n; ++j) {
						LU(k, j) -= LU(k, i) * LU(i, j); //Обновляем матрицу
					}
				}
			}
			return det;
		}
		// Обратная матрица
		Matrix inverse(double tol = 1e-12) const {
			if (rows != cols)
				throw invalid_argument("Matrix must be square to compute inverse");

			int n = rows;
			Matrix inv(n, n);  // Обратная матрица
			Matrix LU(*this);   // Копия исходной матрицы
			vector<int> perm(n); // Вектор перестановок строк
			iota(perm.begin(), perm.end(), 0); // Инициализация: [0, 1, 2, ..., n-1]
			for (int i = 0; i < n; ++i) {
				// Поиск максимального элемента в столбце i
				int maxRow = i;
				for (int k = i + 1; k < n; ++k)
					if (abs(LU(k, i)) > abs(LU(maxRow, i)))
						maxRow = k;
				// Перестановка строк, если необходимо
				if (maxRow != i) {
					swap(perm[i], perm[maxRow]);
					for (int k = 0; k < n; ++k)
						swap(LU(i, k), LU(maxRow, k));
				}
				// Проверка на вырожденность
				if (abs(LU(i, i)) < tol)
					throw invalid_argument("Matrix is singular (pivotal element < tolerance)");
				// Вычисление множителей L и обновление U
				for (int k = i + 1; k < n; ++k) {
					LU(k, i) /= LU(i, i);  // Сохраняем множитель для L
					for (int j = i + 1; j < n; ++j)
						LU(k, j) -= LU(k, i) * LU(i, j);  // Обновляем U
				}
			}
			// Решение LUX = I для каждого столбца обратной матрицы
			for (int col = 0; col < n; ++col) {
				// Прямой ход: решение LY = Pb (где b - столбец единичной матрицы)
				for (int i = 0; i < n; ++i) {
					inv(i, col) = (perm[i] == col) ? 1 : 0;  // Учитываем перестановки
					for (int k = 0; k < i; ++k)
						inv(i, col) -= LU(i, k) * inv(k, col);  // Используем L[i][k]
				}

				// Обратный ход: решение UX = Y
				for (int i = n - 1; i >= 0; --i) {
					for (int k = i + 1; k < n; ++k)
						inv(i, col) -= LU(i, k) * inv(k, col);  // Используем U[i][k]
					inv(i, col) /= LU(i, i);  // Делим на диагональный элемент U
				}
			}

			return inv;
		}
};

// Итерации Арнольда
void arnoldi(const Matrix& A, int m, Matrix& V, Matrix& H) {
	int n = A.getRows();
	V = Matrix(n, m + 1);// Базисные векторы
	H = Matrix(m + 1, m);// Выходная матрица

	V(0, 0)=1; // Единичный первый базисный вектор
	for (int j = 0; j < m; ++j) {
		// Умножение матрицы на вектор (A*v_j)
		Matrix w(n, 1);
		for (int i = 0; i < n; ++i) {
			double sum = 0;
			for (int k = 0; k < n; ++k) sum += A(i, k) * V(k, j);
			w(i, 0) = sum;
		}
		// Вычисление H и w
		for (int i = 0; i <= j; i++) {
			H(i, j) = 0;
			for (int k = 0; k < n; ++k) H(i, j) += w(k, 0) * V(k, i);
			for (int k = 0; k < n; ++k) w(k, 0) -= H(i, j) * V(k, i);
		}
		// Нормализация
		double norm = 0;
		for (int k = 0; k < n; ++k) norm+=w(k,0)*w(k,0);
		norm=sqrt(norm);
		H(j + 1, j) = norm;
		if (H(j + 1, j) < 1e-10) break;
		for (int i = 0; i < n; ++i) V(i, j + 1) = w(i, 0) / norm;
	}
}
// Вычисление сдвигов
Matrix computeADIShifts(const Matrix& A, int k, int m = 20) {
	int n = A.getRows();
	Matrix shifts(k, 1);
	Matrix V_A, H_A, V_Ainv, H_Ainv;
	Matrix A_inv;
	// Если не удалось вычислить обратную, то берем сдвиги по 0.5
	try {
		A_inv = A.inverse();
	} catch (...) {
		for (int i = 0; i < k; ++i) shifts(i, 0) = 0.5;
		return shifts;
	}
	// Итерации Арнольда
	arnoldi(A, m, V_A, H_A);
	arnoldi(A_inv, m, V_Ainv, H_Ainv);

	vector<double> all_shifts;
	for (int j = 0; j < m; ++j) {
		if (abs(H_A(j, j)) > 1e-12 && !isnan(-1.0 / H_A(j, j)))
			all_shifts.push_back(-1.0 / H_A(j, j));
		if (abs(H_Ainv(j, j)) > 1e-12 && !isnan(-1.0 / H_Ainv(j, j)))
			all_shifts.push_back(-1.0 / H_Ainv(j, j));
	}
	// Если сдвигов нет, берем по 0.5
	if (all_shifts.empty()) {
		for (int i = 0; i < k; ++i) shifts(i, 0) = 0.5;
		return shifts;
	}

	sort(all_shifts.begin(), all_shifts.end());
	// Добавляем минимальный и максимальные сдвиги
	int s = 0;
	shifts(s++, 0) = all_shifts.front();
	if (k > 1 && all_shifts.size() > 1) shifts(s++, 0) = all_shifts.back();
	// Добавляем оставшиеся максимально разнесенные сдвиги
	while (s < k) {
		double max_min_dist = -1;
		int best_idx = -1;
		int l=all_shifts.size();
		for (int i = 0; i<l; ++i) {
			double min_dist = numeric_limits<double>::max();
			for (int j = 0; j < s; ++j) //Ищем мин. расстояние до выбранных сдвигов
				min_dist = min(min_dist, abs(all_shifts[i] - shifts(j, 0)));
			// Если элемент далеко от остальных сдвигов, запоминаем его
			if (min_dist > max_min_dist) {
				max_min_dist = min_dist;
				best_idx = i;
			}
		}
		//Лучший сдвиг добавляем
		if (best_idx >= 0) shifts(s++, 0) = all_shifts[best_idx];
		else shifts(s++, 0) = 0.5;
	}
	return shifts;
}
// Решение уравнения Риккати
Matrix solveRiccatiADI(const Matrix& A, const Matrix& B, const Matrix& Q,
                       const Matrix& R, const Matrix& E, int max_iter = 100,
                       double tau = 0.01, double tol = 1e-8) {
	const int n = A.getRows();
	Matrix P(n, n); // Начальное приближение
	double norm;
	// Считаем сразу некоторые матрицы для упрощения
	Matrix E_trans = E.transpose();
	Matrix B_trans = B.transpose();
	Matrix R_inv = R.inverse();
	Matrix shifts = computeADIShifts(A, max_iter);

	for (int k = 0; k < max_iter; ++k) {
		Matrix P_prev = P;

		double alpha = shifts(k, 0);
		double beta = 1.0 - alpha;

		// Первый полушаг
		Matrix M1 = E + alpha * tau * A;
		// Регуляризация
		double epsilon = 1e-8;
		for (int i = 0; i < M1.getRows(); ++i) {
			M1(i, i) += epsilon;
		}
		Matrix M1_inv = M1.inverse();
		Matrix RHS1 = E_trans * P_prev * E - tau * Q + tau * (E_trans * P_prev * B) * R_inv * (B_trans * P_prev * E);
		Matrix P_half = (M1_inv.transpose() * RHS1) * M1_inv;

		// Второй полушаг
		Matrix M2 = E + beta * tau * A;
		for (int i = 0; i < M2.getRows(); ++i) {
			M2(i, i) += epsilon;
		}
		Matrix M2_inv = M2.inverse();
		Matrix RHS2 = E_trans * P_prev * E - tau * Q + tau * (E_trans * P_half * B) * R_inv * (B_trans * P_prev * E);
		P = (M2_inv.transpose() * RHS2) * M2_inv;

		// Проверка сходимости
		norm = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				norm += pow(P(i, j) - P_prev(i, j), 2);
		norm = sqrt(norm);

		if (norm < tol) {
			break;
		}
	}

	cout<<"Convergence: "<<norm<<endl;
	return P;
}

void writeMatrixToFile(const Matrix& matrix, const string& filename) {
	ofstream file(filename);
	if (!file) {
		throw runtime_error("Failed to open file for writing: " + filename);
	}

	int rows = matrix.getRows();
	int cols = matrix.getCols();

	// Сначала определяем максимальную ширину элемента
	int max_width = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			double val = abs(matrix(i, j)) < 1e-10 ? 0 : matrix(i, j);
			ostringstream oss;
			if (val == floor(val)) oss << (int)val;
			else oss << val;
			max_width = max(max_width, (int)oss.str().length());
		}
	}

	// Теперь записываем матрицу в файл
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			double val = abs(matrix(i, j)) < 1e-10 ? 0 : matrix(i, j);
			ostringstream oss;
			if (val == floor(val)) {
				file << setw(max_width + 1) << (int)val;
			} else {
				file << setw(max_width + 1) << val;
			}
		}
		file << endl;
	}

	file.close();
	cout << "Result saved to " << filename << endl;
}
int main() {
	string w="100";
	cout<<w+"x"+w+"\n";
	try {
		Matrix A, B, C, E;
		string u="EABC_tests\\"+w+"x"+w+"\\";
		A.readFromFile(u+"A.dat");
		B.readFromFile(u+"B.dat");
		C.readFromFile(u+"C.dat");
		E.readFromFile(u+"E.dat");
		int t = B.getCols();
		Matrix R(t, t);
		for (int i = 0; i < t; ++i) R(i, i) = 1;
		Matrix Q = C.transpose() * C;

		int n = A.getRows();
		double tau = 0.0001;
		Matrix P = solveRiccatiADI(A, B, Q, R, E, 100, tau);
		cout << "Solution P:\n" << P << endl;
		writeMatrixToFile(P, "P_result_k+1_k_"+w+"x"+w+".dat");

	} catch (const exception& e) {
		cerr << "Error: " << e.what() << endl;
		return 1;
	}

	return 0;
}
