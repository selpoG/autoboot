using System.Collections.Generic;
using System.Text;
using Sprache;

namespace GAPToMathematica
{
	class IrreducibleRepresentations
	{
		public readonly RepresentativesOfConjugacyClasses Representatives;
		public readonly List<List<AlgebraicExpression[,]>> Irreps;
		public IrreducibleRepresentations(RepresentativesOfConjugacyClasses rep, IEnumerable<List<List<List<AlgebraicExpression>>>> irreps)
		{
			Representatives = rep;
			Irreps = new List<List<AlgebraicExpression[,]>>();
			foreach (var ir in irreps)
			{
				var tmp = new List<AlgebraicExpression[,]>();
				foreach (var r in ir)
				{
					var size = r.Count;
					var mat = new AlgebraicExpression[size, size];
					for (var i = 0; i < size; i++) for (var j = 0; j < size; j++) mat[i, j] = r[i][j];
					tmp.Add(mat);
				}
				Irreps.Add(tmp);
			}
		}
		public override string ToString()
		{
			var sb = new StringBuilder("{");
			foreach (var irrep in Irreps)
			{
				sb.Append("{");
				foreach (var mat in irrep)
				{
					var s = mat.GetLength(0);
					sb.Append("{");
					for (var i = 0; i < s; i++)
					{
						sb.Append("{");
						for (var j = 0; j < s; j++)
						{
							if (j != s - 1) sb.Append($"{mat[i, j]}, ");
							else sb.Append(mat[i, j] + "}");
						}
						if (i != s - 1) sb.Append(", ");
						else sb.Append("}");
					}
					sb.Append(", ");
				}
				sb.Remove(sb.Length - 2, 2);
				sb.Append("}, ");
			}
			sb.Remove(sb.Length - 2, 2);
			sb.Append("}");
			return sb.ToString();
		}
		public static readonly Parser<IrreducibleRepresentations> Parser = from rep in RepresentativesOfConjugacyClasses.Parser
																		   from mark in Parse.String("(irrep)=")
																		   from irreps in (from mats in AlgebraicExpression.Parser.ToSequence().ToSequence().ToSequence()
																						   select mats).AtLeastOnce()
																		   select new IrreducibleRepresentations(rep, irreps);
	}
}
