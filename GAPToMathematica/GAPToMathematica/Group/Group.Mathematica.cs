using System.IO;
using System.Linq;
using System.Text;

namespace GAPToMathematica
{
	partial class Group
	{
		public void SaveAsMathematica() => SaveAsMathematica($"sg.{Id.Order}.{Id.Id}.m");
		public void SaveAsMathematica(string path) => File.WriteAllText(path, ToMathematica());
		public string ToMathematica() => string.Format(mathematicaFormat,
													   Id.Order, Id.Id,
													   Sizes,
													   Characters,
													   Generators,
													   Irreps);
		static Group()
		{
			var sb = new StringBuilder();
			for (var i = 0; i < mathematicaFormat.Length; i++)
			{
				var c = mathematicaFormat[i];
				if (c == '$') sb.Append("{" + mathematicaFormat[++i] + "}");
				else
				{
					sb.Append(c);
					if (c == '{' || c == '}') sb.Append(c);
				}
			}
			mathematicaFormat = sb.ToString();
		}
		static readonly string mathematicaFormat = @"data[
{$0, $1},
$2,
$3,
$4,
$5
]
";
	}
}
