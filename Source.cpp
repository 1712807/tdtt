#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <float.h>
#include <climits>
using namespace std;

// bai 1 // nlog(n)
int count(int a[], int x, int n)
{	
	sort(a, a + n);
	int sum = 0;
	int i = 0, j = n - 1;
	while(i < j)
	{
		if (a[i] + a[j] <= x)
		{
			sum += j - i;
			i++;
		}
		else
		{
			j--;
		}
	}
	return sum;
}

// bai 2 // chua toi uu O(n^2)
int getSumIndexK(int a[], int n, int k)
{
	sort(a, a + n);
	vector<int> sum;
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			sum.push_back(a[i] + a[j]);
		}
	}

	sort(std::begin(sum), std::end(sum)); // nlog(n)

	return sum[k - 1];
}

// bai 3 // nlog(n)
vector<vector<int>> find3Index(vector<int> a, int n, int x)
{
	sort(a.begin(), a.end());
	vector<vector<int>> res;
	vector<int> tmp;

	for (int i = 0; i < n - 2; i++)
	{
		int j = i + 1, k = n - 1;
		while (j < k)
		{
			if (a[i] + a[j] + a[k] == x)
			{
				tmp.push_back(i);
				tmp.push_back(j);
				tmp.push_back(k);
				res.push_back(tmp);
				tmp.clear();
				j++;
			}
			else if (a[i] + a[j] + a[k] < x)
			{
				j++;
			}
			else {
				k--;
			}
		}
	}
	return res;
}

// bai 4

// bai 5 
vector<int> prefixSum(int a[], int n)
{
	vector<int> b;
	b.push_back(a[0]);
	for (int i = 1; i < n; i++)
	{
		b.push_back(b[i - 1] + a[i]);
	}
	return b;
}

int CountSum(int a[], int n, int x)
{
	vector<int> prefix = prefixSum(a, n);
	int count = 0;
	for (int i = 0; i < n; i++)
	{
		if (prefix[i] == x) count++;
	}

	sort(std::begin(prefix), std::end(prefix));

	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (prefix[j] - prefix[i] == x)
			{
				count++;
			}
		}
	}
	
	/*int l = 0, r = n - 1;
	while (l < r)
	{
		int i = l, j = r, m = (i + j) / 2;
		
		if (prefix[r] - prefix[m] == x)
		{
			count++;
			prefix.erase(prefix.begin() + m);

			l = 0; r = prefix.size() - 1;
		}
		else if (prefix[r] - prefix[m] <  x)
		{
			r = m;
		}
		else {
			l = m + 1;
		}
	}*/

	return count;

}

// bai 6
int countSubArrDes(int a[], int n)
{
	int *d = new int[n];
	int *s = new int[n];
	for (int i = 0; i < n; i++)
	{
		d[i] = 1;
		s[i] = 0;
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < i; j++)
			if ((a[i] < a[j]) && (d[i] < d[j] + 1))  d[i] = d[j] + 1;

	int max = 0;
	for (int i = 1; i < n; i++)
		if (d[max] < d[i]) max = i;

	for (int i = 0; i < n; i++)
	{
		if (d[i] == 1) s[i] = 1;
		else
			for (int j = 0; j < i; j++)
				if ((a[j] > a[i]) && (d[i] == d[j] + 1)) s[i] = s[i] + s[j];
	}
	return s[max];
}
// bai 7
int findBinary(int a[], int n, int x)
{
	int l = 0, r = n - 1;
	int m = -1;
	while (l < r)
	{
		m = (l + r) / 2;
		if (a[m] == x) return m;
		else if (a[m] < x)
		{
			l = m + 1;
		}
		else r = m;
	}

	return l;
}

int countPair(int a[], int b[], int nA, int nB)
{
	int count = 0;
	sort(a, a + nA);
	sort(b, b + nB);

	for (int i = 0; i < nA; i++)
	{
		if (a[i] < b[0]) count += nB;
		if (a[i] > b[nB - 1]) continue;
		if (a[i] >= b[0] && a[i] <= b[nB - 1]) // tim a[i] trong b[]
		{
			int index = findBinary(b,nB, a[i]);
			if (b[index] != a[i]) count += nB - index;
			else
			{
				for (int j = index + 1; j < nB; j++)
				{
					if (b[j] == a[i]) index = j;
					else break;
				}
				count += nB - 1 - index;
			}
		}
	}


	return count;
}

// bai 8 // ai + bj < s
int countPairSmallerS(int a[], int b[], int nA, int nB, int s)
{
	int count = 0;
	sort(a, a + nA);
	sort(b, b + nB);

	for (int i = 0; i < nA; i++)
	{
		int index = findBinary(b, nB, s - a[i]);
		for (int j = index - 1; j >= 0; j--)
		{
			if (b[j] == s - a[i]) index = j;
			else break;
		}
		cout << index << endl;
		count += index;
	}
	return count;
}

// bai 9 max - ai + bj < s
int maxPairSmallerS(int a[], int b[], int nA, int nB, int s)
{
	int maxSum = -UINT16_MAX;
	sort(a, a + nA);
	sort(b, b + nB);

	for (int i = 0; i < nA; i++)
	{
		int bj = s - a[i];
		int index = findBinary(b, nB, bj);
		
		if (bj > b[nB - 1])
		{
			maxSum = max(maxSum, a[i] + b[nB - 1]);
		}
		else if (bj < b[0])
		{
			continue;
		}
		else
		{
			maxSum = max(maxSum, a[i] + b[index - 1]);
		}
	}
	return maxSum;
}

// bai 10 ai + bj = x
int countPairEqualX(int a[], int b[], int nA, int nB, int x)
{
	int count = 0;
	sort(a, a + nA);
	sort(b, b + nB);

	for (int i = 0; i < nA; i++)
	{
		int bj = x - a[i];
		int index = findBinary(b, nB, x - a[i]);
		
		if (b[index] != bj)
		{
			continue;
		}
		else
		{
			count++;
			for (int j = index + 1; j < nB; j++)
			{
				if (b[j] == bj) count++;
				else break;
			}
			for (int j = index - 1; j > 0; j--)
			{
				if (b[j] == bj) count++;
				else break;
			}
		}
	}

	return count;
}

// bai 11 // nlogn
vector<int> findTheSameValue(int a[], int b[], int nA, int nB)
{
	sort(a, a + nA);
	sort(b, b + nB);

	vector<int> res;
	for (int i = 0; i < nA; i++)
	{
		int index = findBinary(b, nB, a[i]);
		if (b[index] != a[i])
		{
			res.push_back(a[i]);
		}
	}
	return res;
}

// bai 13 nlogn
int getMinSub(int a[], int b[], int nA, int nB)
{
	sort(a, a + nA);
	sort(b, b + nB);
	int minSub = abs(a[0] -b[0]);

	
	for (int i = 0; i < nA; i++)
	{
		int index = findBinary(b, nB, a[i]);
		if (b[index] == a[i]) return 0; ////a[i] co trong b

		if (a[i] < b[0])
		{
			minSub = min(minSub, abs(a[i] - b[0]));
		}
		else if (a[i] > b[nB - 1])
		{
			minSub = min(minSub, abs(a[i] - b[nB - 1]));
		}
		else
		{
			int mintmp = min(abs(b[index - 1] - a[i]), abs(b[index] - a[i]));
			minSub = min(minSub, mintmp);
		}
	}

	return minSub;
}

// bai 14
int getMinSum(int a[], int b[], int nA, int nB)
{
	sort(a, a + nA);
	sort(b, b + nB);

	int minSum = abs(a[0] + b[0]);


	for (int i = 0; i < nA; i++)
	{
		int tmp = -a[i];
		int index = findBinary(b, nB, tmp);
		if (b[index] == tmp) return 0; //// -a[i] co trong b

		if (tmp < b[0])
		{
			minSum = min(minSum, abs(tmp - b[0]));
		}
		else if (tmp > b[nB - 1])
		{
			minSum = min(minSum, abs(tmp - b[nB - 1]));
		}
		else
		{
			int mintmp = min(abs(b[index - 1] - tmp), abs(b[index] - tmp));
			minSum = min(minSum, mintmp);
		}
	}

	return minSum;
}

// bai 15 - O(n)
int findsubArrMax(int a[], int n)
{
	vector<int> s = prefixSum(a, n);

	int minS = s[0];
	int sum_max = s[0];
	int start = 0, end = 0;

	for (int i = 1; i < n; i++)
	{
		if (s[i - 1] < minS)
		{
			minS = s[i - 1];
			start = i - 1;
		}

		if (sum_max < s[i] - minS)
		{
			sum_max = s[i] - minS;
			end = i;
		}
	}

	// sub arr
	for (int i = start + 1; i <= end; i++)
	{
		cout << a[i] << " ";
		cout << endl;
	}
	return sum_max;
}

// bai 16 - sum = 0
struct node {
	int S;
	int index;
};
bool cp(node a, node b) { 
	return a.S < b.S; 
}

void findsubArrLargest(int a[], int n)
{
	node* b = new node[n];
	node tmp;
	int s = 0, len = 0, max_len = 0, l = -1, r = -1;
	for (int i = 0; i < n; i++) {
		s += a[i];
		b[i] = node({ s,i });
	}
	for (int i = n - 1; i >= 0; i--) {
		if (b[i].S == 0)
			if (i > max_len) {
				max_len = i;
				l = -1; r = i;
			}
	}
	sort(b, b + n, cp);
	int minIndex = b[0].index, maxIndex = b[0].index;
	long long int Snow = b[0].S;
	for (int i = 0; i < n; i++) {
		if (Snow == b[i].S) {
			if (b[i].index > maxIndex)
				maxIndex = b[i].index;
			if (b[i].index < minIndex)
				minIndex = b[i].index;
		}
		else {
			len = maxIndex - minIndex;
			if (len > max_len) {
				max_len = len;
				l = minIndex;
				r = maxIndex;
			}
			Snow = b[i].S;
			minIndex = maxIndex = b[i].index;
		}
	}
	len = maxIndex - minIndex;
	if (len > max_len) {
		max_len = len;
		l = minIndex;
		r = maxIndex;
	}
	
	for (int i = l + 1; i <= r; i++)
	{
		cout << a[i] << " ";
	}

	cout << endl;


}

// bai 17 - sum = 0 / nlogn
int countsubArr_Sum0(int a[], int n)
{
	vector<int> s = prefixSum(a, n);
	int count = 0;
	for (int i = 1; i < n; i++)
	{
		if (s[i] == 0) count++;
	}

	sort(s.begin(), s.end());

	int d = 1;
	for (int i = 1; i < n; i++)
	{
		if (s[i] == s[i - 1])
		{
			d++;
		}
		else
		{
			count += d * (d - 1) / 2;
			d = 1;
		}
	}
	count += d * (d - 1) / 2;

	return count;
}

// bai 18 - sum max with min k element
int maxSumK(int a[], int n, int k)
{
	int *maxSum = new int[n];
	maxSum[0] = a[0];

	int curr_max = a[0];
	for (int i = 1; i < n; i++)
	{
		curr_max = max(a[i], curr_max + a[i]);
		maxSum[i] = curr_max;
	}
	
	int sum = 0;
	for (int i = 0; i < k; i++)
		sum += a[i];

	int result = sum;
	for (int i = k; i < n; i++)
	{
		sum = sum + a[i] - a[i - k];
		result = max(result, sum);
		result = max(result, sum + maxSum[i - k]);
	}
	return result;
}


// bai 19
int search(int start, int end, int x, vector<int> p)
{
	if (p[start+1] - p[start] == x)
		return start;

	int l = start;
	int r = end;
	while (r - l > 1)
	{
		int mid = (l + r) / 2;
		int s = p[mid+1] - p[start];
		if (s > x)
		{
			r = mid;
		}
		else if (s < x) l = mid;
		else return mid+1;
	}
	if (p[r+1] - p[start] == x)
		return r;
	else return -1;
}
void findSubArrShortestSumX(int a[], int n, int x)
{
	int *b = new int[2 * n];
	for (int i = 0; i < 2 * n; i++)
	{
		b[i] = a[i%n];
	}
	vector<int> sum = prefixSum(b, 2 * n);

	// o element
	if (x < sum[0] || x > sum[sum.size() - 1]) return;
	// 1 element
	for (int i = 0; i < n; i++)
		if (a[i] == x)
		{
			cout << a[i] << endl;
			return;
		}

	///////
	int s = 0, e = n-1;
	for (int i = 1; i < n; i++)
	{
		int p = search(i, i + n-1, x, sum);
		if (p != -1)
		{
			if (p - i < e - s)
			{
				s = i;
				e = p;
			}
		}

	}

	for (int i = s+1; i <= (e%n); i++)
		cout << a[i] << " ";
}
// bai 20

// bai 21
int maxSumRec(int n, int m)
{
	int a[4][6] = { {8,-9,-4,9,-3,-9},
	{-8,-1,-1,-8,-7,-10},
	{-8,-10,0,-2,-10,6},
	{-3,-5,2,9,-1,1} };

	int prefix_sum[500][500] = {0};
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (j == 0) prefix_sum[i][j] = a[i][j];
			else prefix_sum[i][j] = prefix_sum[i][j - 1] + a[i][j];
		}
	}

	int res = -INFINITY;
	for (int l = 0; l < m; l++) {
		for (int r = l; r < m; r++) {
			int tmp = 0;
			for (int i = 0; i < n; i++) {
				if (l == 0) tmp += prefix_sum[i][r];
				else tmp += prefix_sum[i][r] - prefix_sum[i][l - 1];
				res = max(res, tmp);
				
			}
		}
	}
	return res;
}
// bai 22

class Point_2D
{
public:
	int m_x, m_y, m_idx;

public:
	static Point_2D* O;
	Point_2D(int x, int y) : m_x{ x }, m_y{ y }{}

	bool compareSlope(Point_2D* A, Point_2D* B)
	{
		double slope_a, slope_b;
		if (A->m_x == Point_2D::O->m_x)
		{
			if (A->m_y - Point_2D::O->m_y > 0)
				slope_a = DBL_MAX;
			else
				slope_a = -DBL_MAX;
		}
		else
			slope_a = (A->m_y - Point_2D::O->m_y) * 1.0 / (A->m_x - Point_2D::O->m_x);

		if (B->m_x == Point_2D::O->m_x)
		{
			if (B->m_y - Point_2D::O->m_y > 0)
				slope_b = DBL_MAX;
			else
				slope_b = -DBL_MAX;
		}
		else
			slope_b = (B->m_y - Point_2D::O->m_y) * 1.0 / (B->m_x - Point_2D::O->m_x);

		return slope_a > slope_b;
	}
};

Point_2D* Point_2D::O = new Point_2D(0, 0);

//int main() {
//	/* Enter your code here. Read input from STDIN. Print output to STDOUT */
//
//	int n, x, y, min_x = INT_MAX, pos_O;
//	vector <Point_2D*> map;
//	Point_2D* tree_O = nullptr, *tree_O2;
//
//	n = 6;
//	int a[6] = { 2,9,2,9,11,4 };
//	int b[6] = { 1,6,4,2,4,7 };
//	for (int i = 0; i < n; ++i)
//	{
//		map.push_back(new Point_2D(a[i], b[i]));
//
//		// tim O co hoanh do x nho nhat
//		if (a[i] < min_x)
//		{
//			min_x = a[i];
//			pos_O = i;
//		}
//	}
//
////	 loai diem tree_O ra khoi mang, va tree_O tro thanh goc
//	tree_O = map[pos_O];
//	if (pos_O != n - 1)
//	{
//		map[pos_O] = map[n - 1];
//	}
//	map.pop_back();
//	Point_2D::O = tree_O;
//
//	// sap xep mang theo do doc cua tree_O
//	sort(map.begin(), map.end(), Point_2D::compareSlope);
//
//	// diem trung vi tree_O2
//	tree_O2 = map[map.size() / 2];
//
//	cout << tree_O->m_x << ' '<< tree_O->m_y <<" -> "<< tree_O2->m_x << " "<< tree_O2->m_y << endl;
//
//	for (int i = 0; i < n - 1; ++i)
//	{
//		delete map[i];
//	}
//	delete tree_O;
//	return 0;
//}

// bai 23
struct Point {
	int x;
	int y;
};
int S(Point A, Point B, Point C)
{
	return (B.y + A.y)*(B.x - A.x) + (C.y + B.y)*(C.x - B.x) + (A.y + C.y)*(A.x - C.x);
}
vector<Point> insidePolygon(vector<Point> A, vector<Point> P)
{
	vector<Point> res;
	int n = A.size();
	int m = P.size();
	bool out = false;
	for (int i = 0; i < m; i++)
	{
		int preS = 1;
		bool out = false;

		for (int j = 0; j < n; j++)
		{
			int next = j < n - 1 ? j + 1 : 0;
			int s = S(A[j], A[next], P[i]);

			if (s == 0 || (s*preS < 0 && j>0))
			{
				out = true;
				break;
			}
			preS = s;
		}

		if (!out) res.push_back(P[i]);
	}

	return res;
}

void main()
{
	
	vector<int> a = { 7,5,9,6,2,1,3 };

	vector<vector<int>> b = find3Index(a, a.size(), 10);
	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < 3; j++)
			cout << b[i][j] << " ";

		cout << endl;
	}
	system("pause");
}



