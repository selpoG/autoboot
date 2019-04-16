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
		public static Group CyclicGroup(int n)
		{
			if (n <= 1) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"CyclicGroup({n})"));
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
		public static Group ElementaryAbelianGroup(int n)
		{
			var f = n.Factorize();
			if (f.Count != 1) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"ElementaryAbelianGroup({n})"));
		}
		public static Group DihedralGroup(int n)
		{
			if (n < 2 || n % 2 != 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"DihedralGroup({n})"));
		}
		public static Group QuaternionGroup(int n)
		{
			if (n < 4 || n % 4 != 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"QuaternionGroup({n})"));
		}
		public static Group DicyclicGroup(int n)
		{
			if (n < 4 || n % 4 != 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"DicyclicGroup({n})"));
		}
		public static Group ExtraspecialGroup(int n, int e)
		{
			var f = n.Factorize();
			if (f.Count != 1) throw new ArgumentOutOfRangeException();
			var p = f.Keys.First();
			var q = f[p];
			if (p == 2 || q % 2 == 0) throw new ArgumentOutOfRangeException();
			if (e != p && e != p * p) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"ExtraspecialGroup({n}, {e})"));
		}
		public static Group ExtraspecialGroup(int n, bool plus)
		{
			var f = n.Factorize();
			if (f.Count != 1) throw new ArgumentOutOfRangeException();
			var p = f.Keys.First();
			var q = f[p];
			if (q % 2 == 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"ExtraspecialGroup({n}, {(plus ? "'+'" : "'-'")})"));
		}
		public static Group AlternatingGroup(int n)
		{
			if (n < 3) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"AlternatingGroup({n})"));
		}
		public static Group SymmetricGroup(int n)
		{
			if (n < 2) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"SymmetricGroup({n})"));
		}
		public static Group MathieuGroup(int n)
		{
			if (n < 9 || (12 < n && n < 21) || 24 < n) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"MathieuGroup({n})"));
		}
		public static Group SuzukiGroup(int n)
		{
			var f = n.Factorize();
			if (f.Count != 1 || !f.ContainsKey(2) || f[2] % 2 == 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"SuzukiGroup({n})"));
		}
		public static Group ReeGroup(int n)
		{
			var f = n.Factorize();
			if (f.Count != 1 || !f.ContainsKey(3) || f[3] % 2 == 0) throw new ArgumentOutOfRangeException();
			return GroupFromGAP(Func.GroupInfo($"ReeGroup({n})"));
		}
		public static Group SmallGroup(int order, int id)
		{
			if (order < 2 || id < 1) throw new ArgumentOutOfRangeException();
			var info = Func.GroupInfo(order, id, 300 * 1000);
			if (info == null) throw new TimeoutException();
			var g = GroupFromGAP(info);
			Debug.Assert(g == null || (g.Id.Order == order && g.Id.Id == id));
			return g;
		}
	}
}
