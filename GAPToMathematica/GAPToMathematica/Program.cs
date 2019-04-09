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
			Console.Write("order >> ");
			var G = int.Parse(Console.ReadLine());
			var num = Func.NumberOfGroups(G);
			Console.WriteLine($"There are {num} groups of order {G}.");
			Console.Write($"id ({1} ~ {num}) >> ");
			var i = int.Parse(Console.ReadLine());
			try
			{
				var g = Group.SmallGroup(G, i);
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
