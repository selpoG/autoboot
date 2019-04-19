using System;
using System.Diagnostics;
using System.Linq;
using System.Collections.Generic;
using Sprache;

namespace GAPToMathematica
{
	partial class Group
	{
		public static Group GroupFromGAP(string info)
		{
			var m = Parser.End().TryParse(info.Replace(" ", "").Replace("\n", ""));
			return m.WasSuccessful ? m.Value : null;
		}
		// C_n = (x^n=e)
		public static Group CyclicGroup(int n, bool unitary = false)
		{
			if (n <= 1) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"CyclicGroup({n})", unitary));
		}
		// C_(n1) * ... * C_(nk)
		public static Group AbelianGroup(params int[] ns)
		{
			var pos = new List<int>();
			foreach (var n in ns) if (n > 1) pos.Add(n);
				else if (n < 1) throw new ArgumentOutOfRangeException();
			if (pos.Count < 1) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"AbelianGroup([{string.Join(", ", pos)}])"));
		}
		public static Group ElementaryAbelianGroup(int n, bool unitary = false)
		{
			var f = n.Factorize();
			if (f.Count != 1) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"ElementaryAbelianGroup({n})", unitary));
		}
		public static Group DihedralGroup(int n, bool unitary = false)
		{
			if (n < 2 || n % 2 != 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"DihedralGroup({n})", unitary));
		}
		public static Group QuaternionGroup(int n, bool unitary = false)
		{
			if (n < 4 || n % 4 != 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"QuaternionGroup({n})", unitary));
		}
		public static Group DicyclicGroup(int n, bool unitary = false)
		{
			if (n < 4 || n % 4 != 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"DicyclicGroup({n})", unitary));
		}
		public static Group ExtraspecialGroup(int n, int e, bool unitary = false)
		{
			var f = n.Factorize();
			if (f.Count != 1) throw new ArgumentOutOfRangeException();
			var p = f.Keys.First();
			var q = f[p];
			if (p == 2 || q % 2 == 0) throw new ArgumentOutOfRangeException();
			if (e != p && e != p * p) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"ExtraspecialGroup({n}, {e})", unitary));
		}
		public static Group ExtraspecialGroup(int n, bool plus, bool unitary = false)
		{
			var f = n.Factorize();
			if (f.Count != 1) throw new ArgumentOutOfRangeException();
			var p = f.Keys.First();
			var q = f[p];
			if (q % 2 == 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"ExtraspecialGroup({n}, {(plus ? "'+'" : "'-'")})", unitary));
		}
		public static Group AlternatingGroup(int n, bool unitary = false)
		{
			if (n < 3) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"AlternatingGroup({n})", unitary));
		}
		public static Group SymmetricGroup(int n, bool unitary = false)
		{
			if (n < 2) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"SymmetricGroup({n})", unitary));
		}
		public static Group MathieuGroup(int n, bool unitary = false)
		{
			if (n < 9 || (12 < n && n < 21) || 24 < n) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"MathieuGroup({n})", unitary));
		}
		public static Group SuzukiGroup(int n, bool unitary = false)
		{
			var f = n.Factorize();
			if (f.Count != 1 || !f.ContainsKey(2) || f[2] % 2 == 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"SuzukiGroup({n})", unitary));
		}
		public static Group ReeGroup(int n, bool unitary = false)
		{
			var f = n.Factorize();
			if (f.Count != 1 || !f.ContainsKey(3) || f[3] % 2 == 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"ReeGroup({n})", unitary));
		}
		public static Group SmallGroup(int order, int id, bool unitary = false)
		{
			if (order < 2 || id < 1) throw new ArgumentOutOfRangeException();
			var info = Func.GroupInfo(order, id, 300 * 1000, unitary);
			if (info == null) throw new TimeoutException();
			var g = GroupFromGAP(info);
			Debug.Assert(g == null || (g.Id.Order == order && g.Id.Id == id));
			return g;
		}
	}
}
