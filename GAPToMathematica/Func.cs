using System;
using System.IO;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Sprache;

namespace GAPToMathematica
{
	static class Func
	{
		public static Parser<T> Or<T>() => x => new Result<T>(default(T), x);
		public static Parser<T> Or<T>(this Parser<T> p, params Parser<T>[] qs)
		{
			foreach (var q in qs) p = Parse.Or(p, q);
			return p;
		}
		public static IEnumerable<T> Concat<T>(this T x, IEnumerable<T> y)
		{
			yield return x;
			foreach (var t in y) yield return t;
		}
		public static IEnumerable<T> Concat<T>(this IEnumerable<T> x, T y)
		{
			foreach (var t in x) yield return t;
			yield return y;
		}
		public static Parser<char> Sign = from t in Parse.Char(c => c == '+' || c == '-', "+-")
										  select t;
		public static Parser<char> LowerAlphabet = from t in Parse.Char(c => 'a' <= c && c <= 'z', "a-z")
												   select t;
		public static Parser<char> UpperAlphabet = from t in Parse.Char(c => 'A' <= c && c <= 'Z', "A-Z")
												   select t;
		public static Parser<char> Alphabet = from t in Parse.Char(c => ('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z'), "A-Za-z")
											  select t;
		public static Parser<char> HeadDigit = from t in Parse.Char(c => '1' <= c && c <= '9', "1-9")
											   select t;
		public static Parser<char> Digit = from t in Parse.Char(c => '0' <= c && c <= '9', "0-9")
										   select t;
		public static Parser<int> PositiveInteger = from t in HeadDigit
													from b in Digit.Many()
													select int.Parse(new string(Concat(t, b).ToArray()));
		public static Parser<int> NonzeroInteger = Or(PositiveInteger,
													  from minus in Parse.Char('-')
													  from t in HeadDigit
													  from b in Digit.Many()
													  select -int.Parse(new string(Concat(t, b).ToArray())));
		public static Parser<int> Integer = Or(NonzeroInteger, from zero in Parse.Char('0') select 0);
		// "[ (comma seperated sequence of p) ]"
		// "[  ]", "[ p ]", "[ p, p ]", ...
		public static Parser<List<T>> ToSequence<T>(this Parser<T> p)
		{
			var pc = from x in p
					 from c in Parse.String(",")
					 select x;
			var none = from oc in Parse.String("[]")
					   select new List<T>();
			var many = from o in Parse.String("[")
					   from x in pc.Many()
					   from y in p
					   from c in Parse.String("]")
					   select Concat(x, y).ToList();
			return none.Or(many);
		}
		public static int NumberOfGroups(int order) => int.Parse(Exec(Path.Combine(E.ExeFileDirectory, "numgroup.sh") + $" {order}"));
		public static string GroupInfo(string group, bool unitary = false)
			=> Exec(Path.Combine(E.ExeFileDirectory, "groupinfo.sh") + $" \\\"{group}\\\"" + (unitary ? " :unitary" : ""));
		public static string GroupInfo(int order, int id, bool unitary = false)
			=> Exec(Path.Combine(E.ExeFileDirectory, "smallgroup.sh") + $" {order} {id}" + (unitary ? " :unitary" : ""));
		public static string GroupInfo(int order, int id, int milliseconds, bool unitary = false)
			=> Exec(Path.Combine(E.ExeFileDirectory, "smallgroup.sh") + $" {order} {id}" + (unitary ? " :unitary" : ""), milliseconds);
		public static string Exec(string command)
		{
			var psi = new ProcessStartInfo
			{
				FileName = "/bin/sh",
				UseShellExecute = false,
				RedirectStandardError = true,
				RedirectStandardOutput = true,
				Arguments = $"-c \"{command}\"",
				WorkingDirectory = E.WorkingDirectory
			};
			using (var p = Process.Start(psi))
			{
				p.OutputDataReceived += (s, e) => { };
				p.ErrorDataReceived += (s, e) => { };
				p.BeginOutputReadLine();
				p.BeginErrorReadLine();
				p.WaitForExit();
			}
			return File.ReadAllText(Path.Combine(E.WorkingDirectory, "out.txt"));
		}
		public static string Exec(string command, int milliseconds)
		{
			var psi = new ProcessStartInfo
			{
				FileName = "/bin/sh",
				UseShellExecute = false,
				RedirectStandardError = true,
				RedirectStandardOutput = true,
				Arguments = $"-c \"{command}\"",
				WorkingDirectory = E.WorkingDirectory
			};
			using (var p = Process.Start(psi))
			{
				p.OutputDataReceived += (s, e) => { };
				p.ErrorDataReceived += (s, e) => { };
				p.BeginOutputReadLine();
				p.BeginErrorReadLine();
				p.WaitForExit(milliseconds);
				if (!p.HasExited)
				{
					p.Kill();
					p.Close();
					return null;
				}
			}
			return File.ReadAllText(Path.Combine(E.WorkingDirectory, "out.txt"));
		}
		public static int[] CheckBrackets(this string s)
		{
			var left = new SortedSet<char> { '(', '[', '{', };
			var right = new SortedSet<char> { ')', ']', '}', };
			var N = s.Length;
			var st = new Stack<Tuple<char, int>>();
			var ans = new int[N];
			for (var i = 0; i < N; i++)
			{
				if (left.Contains(s[i])) st.Push(new Tuple<char, int>(s[i], i));
				else if (right.Contains(s[i]))
				{
					var j = st.Pop().Item2;
					ans[i] = j;
					ans[j] = i;
				}
				else ans[i] = -1;
			}
			return ans;
		}
		static readonly List<int> primes;
		static Func()
		{
			const int m = 5000000;
			primes = new List<int>(664579) { 2 };
			var composites = new bool[m];
			composites[0] = false;
			for (var p = 0; p < m; p++)
			{
				if (!composites[p])
				{
					var pnum = 2 * p + 3;
					primes.Add(pnum);
					for (var k = 3 * p + 3; k < m; k += pnum) composites[k] = true;
				}
			}
		}
		public static int NumOfPrimes() => primes.Count;
		public static int NthPrime(int n) => primes[n];
		public static Dictionary<long, int> Factorize(this int n) => Factorize((long)n);
		public static Dictionary<long, int> Factorize(this long n)
		{
			var d = new Dictionary<long, int>();
			foreach (var p in primes)
			{
				if (p * p > n) break;
				if (n % p == 0)
				{
					d.Add(p, 0);
					while (n % p == 0)
					{
						n /= p;
						d[p]++;
					}
				}
			}
			if (n > 1) d.Add(n, 1);
			return d;
		}
	}
}
