using System;
using System.Diagnostics;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Sprache;

namespace GAPToMathematica
{
	class E
	{
		static void Main()
		{
			var nsb = new StringBuilder("Errors at :\n");
			var tsb = new StringBuilder("Timeout at :\n");
			var ng = new List<GroupID>();
			var tmout = new List<GroupID>();
			for (var G = 4; G <= 6; G++)
			{
				var num = Func.NumberOfGroups(G);
				Console.WriteLine("{0} => {1}", G, num);
				for (var i = 1; i <= num; i++)
				{
					try
					{
						var g = Group.SmallGroup(G, i);
						if (g != null && g.Check())
						{
							Console.WriteLine("[{0}, {1}]", G, i);
							g.SaveAsMathematica();
						}
						else
						{
							Console.WriteLine("Error! => {0}, {1}", G, i);
							ng.Add(new GroupID(G, i));
							nsb.Append((nsb.Length <= 12 ? "" : ", ") + $"[{G}, {i}]");
						}
					}
					catch (TimeoutException)
					{
						Console.WriteLine("Timeout Error! => {0}, {1}", G, i);
						tmout.Add(new GroupID(G, i));
						tsb.Append((tsb.Length <= 13 ? "" : ", ") + $"[{G}, {i}]");
					}
					using (var sw = new StreamWriter("ng.txt"))
					{
						sw.WriteLine(nsb);
						sw.WriteLine(tsb);
					}
				}
			}
			if (ng.Count > 0) Console.WriteLine(nsb);
			if (tmout.Count > 0) Console.WriteLine(tsb);
		}
	}
}
