using System.IO;
using System.Linq;

namespace GAPToMathematica
{
	partial class Group
	{
		public void SaveAsMathematica() => SaveAsMathematica(Path.Combine(E.WorkingDirectory, $"sg.{Id.Order}.{Id.Id}.m"));
		public void SaveAsMathematica(string path) => File.WriteAllText(path, ToMathematica());
		public string ToMathematica() => string.Format(mathematicaFormat,
													   Id,
													   Sizes,
													   Characters,
													   Generators,
													   Irreps,
													   Elements == null ? ""
													   : ",\nFunction[" + Generators.ToString().Replace("\"", "")
													   + ", {" + string.Join(", ", Elements.Select(g => g.ToString(Generators.Generators[0]))).Replace("\"", "") + "}]");
		static readonly string mathematicaFormat = @"data[
{0},
{1},
{2},
{3},
{4}{5}
]
";
	}
}
