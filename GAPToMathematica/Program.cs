using System;
using System.Diagnostics;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Sprache;
using System.Reflection;

namespace GAPToMathematica
{
	class E
	{
		public static readonly string WorkingDirectory = Directory.GetCurrentDirectory();
		public static readonly string ExeFileDirectory = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);
		static void Main()
		{
			Console.Write("unitary? [Y/n] >> ");
			var flag_str = Console.ReadLine();
			var unitary = flag_str.Length == 0 || flag_str[0] == 'Y' || flag_str[0] == 'y';
			Console.WriteLine(unitary ? "unitary irreps" : "non-unitary irreps");
			Console.Write("order >> ");
			var G = int.Parse(Console.ReadLine());
			var num = Func.NumberOfGroups(G);
			Console.WriteLine($"There are {num} groups of order {G}.");
			Console.Write($"id ({1} ~ {num}) >> ");
			var i = int.Parse(Console.ReadLine());
			try
			{
				var g = Group.SmallGroup(G, i, unitary);
				if (g != null && g.Check())
				{
					Console.WriteLine("Success: SmallGroup({0}, {1})", G, i);
					g.SaveAsMathematica();
				}
				else Console.WriteLine($"Error! To debug, please run \"{Path.Combine(ExeFileDirectory, "smallgroup.sh")} {G} {i}\" directly.");
			}
			catch (TimeoutException) { Console.WriteLine("Timeout Error!"); }
		}
	}
}
