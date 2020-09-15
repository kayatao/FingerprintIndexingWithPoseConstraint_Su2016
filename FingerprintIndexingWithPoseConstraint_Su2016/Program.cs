using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.IO;
using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.CV.UI;

namespace FingerprintIndexingWithPoseConstraint_Su2016
{
    class Program
    {
        static void Main(string[] args)
        {
            string projectDirectory = Directory.GetParent(Environment.CurrentDirectory).Parent.Parent.FullName;

            FPIndexingMCCPoseConst_Su2016 su2016 = new FPIndexingMCCPoseConst_Su2016();
            //FPIndexingMCCPoseConst_Su2016 su2016 = new FPIndexingMCCPoseConst_Su2016(projectDirectory + @"\Dataset\Hash.txt");

            //su2016.writeHashFunctionsToText(projectDirectory + @"\Dataset\Hash.txt");


            su2016.enrollFPFromMCCText_AppendixB_WithPoseText(projectDirectory + @"\Dataset\MCCTemplateAppendixB\Gallery",
                                                              projectDirectory + @"\Dataset\Pose");

            //su2016.WriteIndexTableToBinary(projectDirectory + @"\Dataset\IndexTable.bin");
            //su2016.ReadIndexTableFromBinary(projectDirectory + @"\Dataset\IndexTable.bin");

            su2016.queryFPInIndexTableFromMCCText_AppendixB_WithPoseText(projectDirectory + @"\Dataset\MCCTemplateAppendixB\Query",
                                                                         projectDirectory + @"\Dataset\Pose",
                                                                         projectDirectory + @"\Dataset\Result\result.txt",
                                                                         'F', 'S', true,
                                                                         projectDirectory + @"\Dataset\Result", ".png",
                                                                         projectDirectory + @"\Dataset\MCCTemplateAppendixB\Gallery",
                                                                         projectDirectory + @"\Dataset\Pose",
                                                                         projectDirectory + @"\Dataset\FPImage", projectDirectory + @"\Dataset\FPImage", ".bmp");

        }
    }
}
