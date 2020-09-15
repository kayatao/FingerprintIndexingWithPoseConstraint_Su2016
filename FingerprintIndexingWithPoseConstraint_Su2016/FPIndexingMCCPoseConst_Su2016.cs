using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Drawing;
using System.IO;
using System.Collections;
using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.CV.CvEnum;
using Emgu.CV.UI;
using System.Windows.Forms;
using System.Diagnostics;

namespace FingerprintIndexingWithPoseConstraint_Su2016
{
    class FPIndexingMCCPoseConst_Su2016
    {
        #region Data Model
        public class IndexItem
        {
            public string F_ID;
            public int M_ID;

            public double distanceFromPose;
            public double angleFromPose;
            public PointF MInPoseCoordinate;
            public PointF MInFixedPoseCoordinate;
            public double alignedMinutiaDirection;
        }

        public class RetrieveItem : IndexItem
        {
            public int Hits;
        }

        public class Cylinder
        {
            public int M_ID;
            public BitArray cmBit;
            public bool valid = true;
            public PointF M;
            public double direction;

            public double distanceFromPose;
            public double angleFromPose;
            public PointF MInPoseCoordinate;
            public PointF MInFixedPoseCoordinate;
            public double alignedMinutiaDirection;
        };

        public class CylindersTemplate
        {
            public Cylinder[] V;
            public int numOfValidCylinder;
            public Pose pose;
            public Pose alignedPose;
            public double rotateDirectionRadius;
            public double rotateDirectionDegree;
            public string F_ID;
        }

        public class Pose
        {
            public PointF p;
            public double angle;
        }
        #endregion

        #region parameters
        public int Ns = 8;
        public int Nd = 6;
        public int minPC = 2;
        public int p = 30;
        public int h = 32;
        public int hashBits = 4;

        public double eL = 80.0;
        public double eTheta = (20.0 * Math.PI) / 180.0;
        #endregion

        Random rnd;

        List<int[][]> H;
        Dictionary<int, List<IndexItem>>[] indexTable;
        Dictionary<string, int> numOfMForFID;
        Dictionary<string, int> allRankDb;

        public FPIndexingMCCPoseConst_Su2016()
        {
            rnd = new Random();
            indexTable = new Dictionary<int, List<IndexItem>>[h * Nd];
            numOfMForFID = new Dictionary<string, int>();
            allRankDb = new Dictionary<string, int>();

            generateAllHashFunctions();
            initialIndexItemDict();
        }

        public FPIndexingMCCPoseConst_Su2016(string hashTextFilePath, double eL = 80, double eTheta = (20.0 * Math.PI) / 180.0)
        {

            rnd = new Random();
            indexTable = new Dictionary<int, List<IndexItem>>[h * Nd];
            numOfMForFID = new Dictionary<string, int>();
            allRankDb = new Dictionary<string, int>();

            readAllHashFunctionsFromText(hashTextFilePath);
            initialIndexItemDict();
        }

        ~FPIndexingMCCPoseConst_Su2016()
        {
            rnd = null;
            H.Clear();
            H.TrimExcess();
            H = null;
            foreach (Dictionary<int, List<IndexItem>> iTable in indexTable)
            {
                iTable.Clear();
            }
            indexTable = null;
            numOfMForFID.Clear();
            numOfMForFID = null;
            allRankDb.Clear();
            allRankDb = null;
        }

        #region Initial
        private void initialIndexItemDict()
        {
            for (int hf = 0; hf < H.Count; hf++)
            {
                for (int k = 0; k < Nd; k++)
                {
                    indexTable[hf * Nd + k] = new Dictionary<int, List<IndexItem>>();
                }
            }
        }
        #endregion

        #region Get Features From File
        public Pose getPoseFromFile(string path)
        {
            Pose p = new Pose();
            string[] strPoseLines = System.IO.File.ReadAllLines(path);
            char[] splitChars = { '\t', ' ' };
            string[] stringPose = strPoseLines[0].Split(splitChars);
            p.p.X = Convert.ToInt32(stringPose[0]);
            p.p.Y = Convert.ToInt32(stringPose[1]);
            p.angle = (-Math.PI / 2.0) - (Convert.ToDouble(stringPose[2]) * Math.PI / 180.0);

            return p;
        }
        
        public void createBitCylinderTemplateFromMCCFile_AppendixB(string pathFile, ref CylindersTemplate CS, Pose pose)
        {
            string[] strCylinderLines = System.IO.File.ReadAllLines(pathFile);

            Size imgSize = new Size(Convert.ToInt32(strCylinderLines[0]), Convert.ToInt32(strCylinderLines[1]));

            double setPoseDirection = -(Math.PI / 2.0);
            Pose fixedPose = new Pose();
            fixedPose.p = new PointF((float)(imgSize.Width / 2.0), (float)(imgSize.Height / 2.0));
            fixedPose.angle = setPoseDirection;

            double rotatedDirectionRadius = differenceAngleOf2Directions_0To2Pi(fixedPose.angle, pose.angle);
            double rotatedDirectionDegree = rotatedDirectionRadius * 180 / Math.PI;
            double angleP = 0.0;
            double distanceP = 0.0;
            double delta_xP, delta_yP;
            delta_xP = pose.p.X - (imgSize.Width / 2);
            delta_yP = pose.p.Y - (imgSize.Height / 2);
            angleP = Math.Atan2(delta_yP, delta_xP);
            distanceP = Math.Sqrt(Math.Pow(delta_xP, 2.0) + Math.Pow(delta_yP, 2.0));
            double alignedPoseAngle = differenceAngleOf2Directions_0To2Pi(angleP, -rotatedDirectionRadius);
            double alignedPoseDirection = differenceAngleOf2Directions_0To2Pi(pose.angle, -rotatedDirectionRadius);//-(Math.PI);
            float alignedX = (float)((imgSize.Width / 2) + distanceP * Math.Cos(alignedPoseAngle));
            float alignedY = (float)((imgSize.Height / 2) + distanceP * Math.Sin(alignedPoseAngle));
            Pose alignedPose = new Pose();
            alignedPose.p = new PointF(alignedX, alignedY);
            alignedPose.angle = alignedPoseDirection;
            CS.alignedPose = alignedPose;
            CS.rotateDirectionRadius = rotatedDirectionRadius;
            CS.rotateDirectionDegree = rotatedDirectionDegree;

            #region Plot pose on image which is rotated at center of image
            //Image<Bgr, Single> plotMinutiaeWithPose = new Image<Bgr, float>(imgSize.Width, imgSize.Height, new Bgr(Color.White));
            //Image<Bgr, Single> plotMinutiaeWithPoseAlign = plotMinutiaeWithPose.Rotate(rotatedDirectionDegree, new Bgr(Color.Black)).Clone();
            //plotPoseOnImage(pose, ref plotMinutiaeWithPose);
            //plotPoseOnImage(alignedPose, ref plotMinutiaeWithPoseAlign);
            //ImageViewer.Show(plotMinutiaeWithPoseAlign);
            //ImageViewer.Show(plotMinutiaeWithPose.ConcateHorizontal(plotMinutiaeWithPoseAlign));
            #endregion

            int n = 0;
            if (Int32.TryParse(strCylinderLines[3], out n) && n != 0)
            {
                Cylinder[] V = new Cylinder[n];
                int vm = 0;
                for (int ms = 4; ms < 4 + n; ms++)
                {
                    string[] stringMinutiaLocDir = strCylinderLines[ms].Split(' ');
                    V[vm] = new Cylinder();
                    V[vm].M_ID = vm;
                    V[vm].M = new PointF();
                    V[vm].M.X = Convert.ToInt32(stringMinutiaLocDir[0]);
                    V[vm].M.Y = Convert.ToInt32(stringMinutiaLocDir[1]);
                    //angle in MCC format is increasing counter-clockwise system
                    //minus 2Pi to convert angle back in increase clockwise system
                    V[vm].direction = (2 * Math.PI) - Convert.ToDouble(stringMinutiaLocDir[2]);

                    #region plot original minutiae
                    //PointF str1 = new PointF(V[vm].M.X, V[vm].M.Y);
                    //PointF end1 = new PointF();
                    //float length = 15;
                    //end1.X = (float)(str1.X + length * Math.Cos(V[vm].direction));
                    //end1.Y = (float)(str1.Y + length * Math.Sin(V[vm].direction));
                    //LineSegment2DF minutiaLine = new LineSegment2DF(str1, end1);
                    //plotMinutiaeWithPose.Draw(new CircleF(V[vm].M, 5), new Bgr(0, 0, 255));
                    //plotMinutiaeWithPose.Draw(minutiaLine, new Bgr(0, 0, 255), 1);
                    #endregion

                    double angle = 0.0;
                    double distance = 0.0;
                    angle = angleOf2Points(V[vm].M, pose.p);
                    distance = distanceOf2Points(V[vm].M, pose.p);
                    double alignedMinutiaLocationAngle = differenceAngleOf2Directions_0To2Pi(angle, -rotatedDirectionRadius);
                    double alignedMinutiaDirection = differenceAngleOf2Directions_0To2Pi(V[vm].direction, -rotatedDirectionRadius);
                    float alignedMX = (float)(alignedPose.p.X + distance * Math.Cos(alignedMinutiaLocationAngle));
                    float alignedMY = (float)(alignedPose.p.Y + distance * Math.Sin(alignedMinutiaLocationAngle));
                    float alignedFixedMX = (float)(fixedPose.p.X + distance * Math.Cos(alignedMinutiaLocationAngle));
                    float alignedFixedMY = (float)(fixedPose.p.Y + distance * Math.Sin(alignedMinutiaLocationAngle));

                    #region plot aligned minutiae
                    //PointF str2 = new PointF(alignedMX, alignedMY);
                    //PointF end2 = new PointF();
                    //end2.X = (float)(str2.X + length * Math.Cos(alignedMinutiaDirection));
                    //end2.Y = (float)(str2.Y + length * Math.Sin(alignedMinutiaDirection));
                    //LineSegment2DF alignedMinutiaLine = new LineSegment2DF(str2, end2);
                    //plotMinutiaeWithPoseAlign.Draw(new CircleF(str2, 5), new Bgr(0, 0, 255));
                    //plotMinutiaeWithPoseAlign.Draw(alignedMinutiaLine, new Bgr(0, 0, 255), 1);
                    #endregion

                    V[vm].distanceFromPose = distance;
                    V[vm].angleFromPose = angle;
                    V[vm].alignedMinutiaDirection = alignedMinutiaDirection;
                    V[vm].MInPoseCoordinate = new PointF(alignedMX, alignedMY);
                    V[vm].MInFixedPoseCoordinate = new PointF(alignedFixedMX, alignedFixedMY);
                    vm++;
                }
                //ImageViewer.Show(plotMinutiaeWithPose.ConcateHorizontal(plotMinutiaeWithPoseAlign));
                //plotMinutiaeWithPose.Dispose();
                //plotMinutiaeWithPoseAlign.Dispose();

                vm = 0;
                int numOfValidCylinder = 0;
                for (int cm = 4 + n + 1; cm < strCylinderLines.Length; cm++)
                {

                    if (strCylinderLines[cm].Contains("True"))
                    {
                        numOfValidCylinder++;
                        string stringCylinderLine = strCylinderLines[cm].Substring(5); //cut off 'True' word
                        string[] stringCylinderMinutia = stringCylinderLine.Split(' ');

                        BitArray cylinder = new BitArray(Ns * Ns * Nd, false);
                        //Image<Gray, Single>[] imgCm = new Image<Gray, float>[Nd];
                        //int[] cylinderMask = createMaskForHashFunction();

                        for (int k = 0; k < Nd; k++)
                        {
                            //imgCm[k] = new Image<Gray, float>(Ns, Ns);
                            //imgCm[k].SetValue(0.0);
                            for (int i = 0; i < Ns; i++)
                            {
                                for (int j = 0; j < Ns; j++)
                                {
                                    int linString = calLin(i + 1, j + 1, k + 1, Ns) - 1 + (Ns * Ns);
                                    int linCylinder = calLin(i + 1, j + 1, k + 1, Ns) - 1;
                                    cylinder[linCylinder] = stringCylinderMinutia[linString] == "1" ? true : false;

                                    //int linCylinderMask = calLin(i + 1, j + 1, 0 + 1, Ns) - 1;
                                    //double gray = 128;
                                    //if (cylinderMask.Contains(linCylinderMask))
                                    //{
                                    //    gray = cylinder[linCylinder] ? 255.0 : 0.0;
                                    //}
                                    //imgCm[k].Data[j, i, 0] = (float)gray;

                                }
                            }
                        }
                        V[vm].valid = true;
                        V[vm].cmBit = cylinder;

                        //Image<Gray, Single> imgCmAll = imgCm[0].Clone();
                        //imgCm[0].Dispose();
                        //for (int k = 1; k < Nd; k++)
                        //{
                        //    imgCmAll = imgCmAll.ConcateVertical(imgCm[k]);
                        //    imgCm[k].Dispose();
                        //}
                        //ImageViewer.Show(imgCmAll);
                        //imgCmAll.Save("IndexCylinder.png");
                        //imgCmAll.Dispose();
                    }
                    else if (strCylinderLines[cm].Contains("False"))
                    {
                        V[vm].valid = false;
                    }
                    vm++;
                }

                CS.V = (Cylinder[])V.Clone();
                CS.numOfValidCylinder = numOfValidCylinder;
            }
        }
        #endregion

        #region Hash Functions
        private int[] createMaskForHashFunction(double strAngle = -Math.PI, double endAngle = Math.PI)
        {
            //Create Mask
            BitArray cylinderMask = new BitArray(Ns * Ns, false);
            List<int> cellID = new List<int>();
            //Image<Gray, Single> cMask = new Image<Gray, Single>(Ns, Ns);
            //cMask.SetZero();
            for (int i = 0; i < Ns; i++)
            {
                for (int j = 0; j < Ns; j++)
                {
                    int linMask = calLin(i + 1, j + 1, 0 + 1, Ns) - 1;

                    PointF center = new PointF((float)((Ns - 1) / 2.0), (float)((Ns - 1) / 2.0));
                    double distance = distanceOf2Points(new PointF(i, j), center);
                    //double angle = angleOf2Points(center, new PointF(i, j));
                    double angle = angleOf2Points(new PointF(i, j), center);

                    if (angle <= endAngle && angle >= strAngle && distance <= (Ns / 2.0))
                    {
                        cylinderMask[linMask] = true;
                        cellID.Add(linMask);
                    }

                    //cylinderMask[linMask] = (angle <= endAngle && angle >= strAngle) ? true : false;
                    //cylinderMask[linMask] = distance <= (Ns / 2.0) ? true : false;

                    //cMask.Data[j, i, 0] = cylinderMask[linMask] ? (float)255.0 : (float)Double.NaN;

                }
            }
            //ImageViewer.Show(cMask);
            //cMask.Dispose();
            cellID.Sort();
            return cellID.ToArray();
        }

        private int[] generateHashFunctionFromMask(int[] mask, int k)
        {
            int[] hashFunctions = new int[hashBits];

            for (int i = 0; i < hashBits; i++)
            {
                int r = rnd.Next(mask.Length);
                hashFunctions[i] = mask[r] + (k * Ns * Ns);
                mask = mask.Where(m => m != mask[r]).ToArray();
            }

            Array.Sort(hashFunctions);
            return hashFunctions;
        }

        private void generateAllHashFunctions()
        {
            #region create all masks
            int[] maskTLQuarter = createMaskForHashFunction(-Math.PI, -Math.PI / 2.0);
            int[] maskTRQuarter = createMaskForHashFunction(-Math.PI / 2.0, 0);
            int[] maskBLQuarter = createMaskForHashFunction(Math.PI / 2.0, Math.PI);
            int[] maskBRQuarter = createMaskForHashFunction(0, Math.PI / 2.0);

            int[] maskTHalf = new int[maskTLQuarter.Length + maskTLQuarter.Length];
            Array.Copy(maskTLQuarter, maskTHalf, maskTLQuarter.Length);
            Array.Copy(maskTRQuarter, 0, maskTHalf, maskTLQuarter.Length, maskTRQuarter.Length);
            Array.Sort(maskTHalf);

            int[] maskBHalf = new int[maskBLQuarter.Length + maskBRQuarter.Length];
            Array.Copy(maskBLQuarter, maskBHalf, maskBLQuarter.Length);
            Array.Copy(maskBRQuarter, 0, maskBHalf, maskBLQuarter.Length, maskBRQuarter.Length);
            Array.Sort(maskBHalf);

            int[] maskLHalf = new int[maskTLQuarter.Length + maskBLQuarter.Length];
            Array.Copy(maskTLQuarter, maskLHalf, maskTLQuarter.Length);
            Array.Copy(maskBLQuarter, 0, maskLHalf, maskTLQuarter.Length, maskBLQuarter.Length);
            Array.Sort(maskLHalf);

            int[] maskRHalf = new int[maskTRQuarter.Length + maskBRQuarter.Length];
            Array.Copy(maskTRQuarter, maskRHalf, maskTRQuarter.Length);
            Array.Copy(maskBRQuarter, 0, maskRHalf, maskTRQuarter.Length, maskBRQuarter.Length);
            Array.Sort(maskRHalf);

            int[] maskFull = createMaskForHashFunction();
            #endregion

            int[][] H1_4_TL1 = new int[Nd][];
            int[][] H1_4_TL2 = new int[Nd][];

            int[][] H1_4_TR1 = new int[Nd][];
            int[][] H1_4_TR2 = new int[Nd][];

            int[][] H1_4_BL1 = new int[Nd][];
            int[][] H1_4_BL2 = new int[Nd][];

            int[][] H1_4_BR1 = new int[Nd][];
            int[][] H1_4_BR2 = new int[Nd][];

            int[][] H1_2_T1 = new int[Nd][];
            int[][] H1_2_T2 = new int[Nd][];
            int[][] H1_2_T3 = new int[Nd][];

            int[][] H1_2_B1 = new int[Nd][];
            int[][] H1_2_B2 = new int[Nd][];
            int[][] H1_2_B3 = new int[Nd][];

            int[][] H1_2_L1 = new int[Nd][];
            int[][] H1_2_L2 = new int[Nd][];
            int[][] H1_2_L3 = new int[Nd][];

            int[][] H1_2_R1 = new int[Nd][];
            int[][] H1_2_R2 = new int[Nd][];
            int[][] H1_2_R3 = new int[Nd][];

            int[][] HF1 = new int[Nd][];
            int[][] HF2 = new int[Nd][];
            int[][] HF3 = new int[Nd][];
            int[][] HF4 = new int[Nd][];
            int[][] HF5 = new int[Nd][];
            int[][] HF6 = new int[Nd][];
            int[][] HF7 = new int[Nd][];
            int[][] HF8 = new int[Nd][];
            int[][] HF9 = new int[Nd][];
            int[][] HF10 = new int[Nd][];
            int[][] HF11 = new int[Nd][];
            int[][] HF12 = new int[Nd][];

            for (int k = 0; k < Nd; k++)
            {
                H1_4_TL1[k] = generateHashFunctionFromMask(maskTLQuarter, k);
                H1_4_TL2[k] = generateHashFunctionFromMask(maskTLQuarter, k);

                H1_4_TR1[k] = generateHashFunctionFromMask(maskTRQuarter, k);
                H1_4_TR2[k] = generateHashFunctionFromMask(maskTRQuarter, k);

                H1_4_BL1[k] = generateHashFunctionFromMask(maskBLQuarter, k);
                H1_4_BL2[k] = generateHashFunctionFromMask(maskBLQuarter, k);

                H1_4_BR1[k] = generateHashFunctionFromMask(maskBRQuarter, k);
                H1_4_BR2[k] = generateHashFunctionFromMask(maskBRQuarter, k);

                H1_2_T1[k] = generateHashFunctionFromMask(maskTHalf, k);
                H1_2_T2[k] = generateHashFunctionFromMask(maskTHalf, k);
                H1_2_T3[k] = generateHashFunctionFromMask(maskTHalf, k);

                H1_2_B1[k] = generateHashFunctionFromMask(maskBHalf, k);
                H1_2_B2[k] = generateHashFunctionFromMask(maskBHalf, k);
                H1_2_B3[k] = generateHashFunctionFromMask(maskBHalf, k);

                H1_2_L1[k] = generateHashFunctionFromMask(maskLHalf, k);
                H1_2_L2[k] = generateHashFunctionFromMask(maskLHalf, k);
                H1_2_L3[k] = generateHashFunctionFromMask(maskLHalf, k);

                H1_2_R1[k] = generateHashFunctionFromMask(maskRHalf, k);
                H1_2_R2[k] = generateHashFunctionFromMask(maskRHalf, k);
                H1_2_R3[k] = generateHashFunctionFromMask(maskRHalf, k);

                HF1[k] = generateHashFunctionFromMask(maskFull, k);
                HF2[k] = generateHashFunctionFromMask(maskFull, k);
                HF3[k] = generateHashFunctionFromMask(maskFull, k);
                HF4[k] = generateHashFunctionFromMask(maskFull, k);
                HF5[k] = generateHashFunctionFromMask(maskFull, k);
                HF6[k] = generateHashFunctionFromMask(maskFull, k);
                HF7[k] = generateHashFunctionFromMask(maskFull, k);
                HF8[k] = generateHashFunctionFromMask(maskFull, k);
                HF9[k] = generateHashFunctionFromMask(maskFull, k);
                HF10[k] = generateHashFunctionFromMask(maskFull, k);
                HF11[k] = generateHashFunctionFromMask(maskFull, k);
                HF12[k] = generateHashFunctionFromMask(maskFull, k);
            }

            H = new List<int[][]>();

            H.Add(H1_4_TL1);
            H.Add(H1_4_TL2);

            H.Add(H1_4_TR1);
            H.Add(H1_4_TR2);

            H.Add(H1_4_BL1);
            H.Add(H1_4_BL2);

            H.Add(H1_4_BR1);
            H.Add(H1_4_BR2);

            H.Add(H1_2_T1);
            H.Add(H1_2_T2);
            H.Add(H1_2_T3);

            H.Add(H1_2_B1);
            H.Add(H1_2_B2);
            H.Add(H1_2_B3);

            H.Add(H1_2_L1);
            H.Add(H1_2_L2);
            H.Add(H1_2_L3);

            H.Add(H1_2_R1);
            H.Add(H1_2_R2);
            H.Add(H1_2_R3);

            H.Add(HF1);
            H.Add(HF2);
            H.Add(HF3);
            H.Add(HF4);
            H.Add(HF5);
            H.Add(HF6);
            H.Add(HF7);
            H.Add(HF8);
            H.Add(HF9);
            H.Add(HF10);
            H.Add(HF11);
            H.Add(HF12);
        }

        public void writeHashFunctionsToText(string saveFilePath)
        {
            using (var writer = new StreamWriter(saveFilePath))
            {
                int k = 0;
                foreach (int[][] hk in H)
                {
                    writer.WriteLine("{0}", k);
                    foreach (int[] h in hk)
                    {
                        writer.WriteLine("{0}\t{1}\t{2}\t{3}", h[0], h[1], h[2], h[3]);
                    }
                    k++;
                }
            }
        }

        private void readAllHashFunctionsFromText(string saveFilePath)
        {
            H = new List<int[][]>();
            string[] strHLines = System.IO.File.ReadAllLines(saveFilePath);

            for (int l = 1; l < strHLines.Length; l += Nd + 1)
            {
                int[][] tmpH = new int[Nd][];
                for (int n = 0; n < Nd; n++)
                {
                    string[] hk = strHLines[l + n].Split('\t');

                    int[] hashFunctions = new int[hashBits];

                    for (int k = 0; k < hk.Length; k++)
                    {
                        hashFunctions[k] = Convert.ToInt32(hk[k]);
                    }
                    tmpH[n] = hashFunctions;
                }
                H.Add(tmpH);
            }


        }
        #endregion

        #region Enroll and Retrieve
        public void enrollFPFromMCCText_AppendixB_WithPoseText(string galleryMCCFolder, string galleryPoseFolder)
        {
            var vEnrollMCCFileNames = new DirectoryInfo(galleryMCCFolder).GetFileSystemInfos("*.txt").Select(x => x.FullName);
            string[] enrollMCCFileNames = vEnrollMCCFileNames.ToArray();

            string MCCFpath = Path.GetDirectoryName(enrollMCCFileNames[0]);

            var vEnrollPoseFileNames = new DirectoryInfo(galleryPoseFolder).GetFileSystemInfos("*.txt").Select(x => x.FullName);
            string[] enrollPoseFileNames = vEnrollPoseFileNames.ToArray();

            string PoseFpath = Path.GetDirectoryName(enrollPoseFileNames[0]);

            Stopwatch stopwatch = new Stopwatch();

            for (int f = 0; f < enrollMCCFileNames.Length; f++)//107
            {
                // Begin timing.
                stopwatch.Restart();

                string fname = Path.GetFileNameWithoutExtension(enrollMCCFileNames[f]);

                Pose pose = getPoseFromFile(PoseFpath + "//" + fname + ".txt");

                CylindersTemplate CS = new CylindersTemplate();
                CS.F_ID = fname;
                CS.pose = pose;
                createBitCylinderTemplateFromMCCFile_AppendixB(enrollMCCFileNames[f], ref CS, pose);

                enrollToIndexTable(CS);
                numOfMForFID.Add(CS.F_ID, CS.numOfValidCylinder);

                // Stop timing.
                stopwatch.Stop();

                Console.WriteLine("{0} is enrolled successfully. Time elapsed: {1}", fname, stopwatch.Elapsed);

            }
        }

        private void enrollToIndexTable(CylindersTemplate CS)
        {
            for (int vm = 0; vm < CS.V.Length; vm++)
            {
                if (CS.V[vm].valid)
                {
                    int idHashTable = 0;
                    foreach (int[][] hf in H)
                    {
                        for (int k = 0; k < Nd; k++)
                        {
                            BitArray indexTermBits = new BitArray(hashBits);

                            for (int hb = 0; hb < hashBits; hb++)
                            {
                                int getElementAt = hf[k][hb];//H1_4_TL1[k][h];
                                indexTermBits[hb] = CS.V[vm].cmBit[getElementAt];
                            }

                            if (popCount(indexTermBits) >= minPC)
                            {
                                IndexItem idxItem = new IndexItem();
                                idxItem.F_ID = CS.F_ID;
                                idxItem.M_ID = CS.V[vm].M_ID;
                                idxItem.distanceFromPose = CS.V[vm].distanceFromPose;
                                idxItem.angleFromPose = CS.V[vm].angleFromPose;
                                idxItem.MInPoseCoordinate = CS.V[vm].MInPoseCoordinate;
                                idxItem.MInFixedPoseCoordinate = CS.V[vm].MInFixedPoseCoordinate;
                                idxItem.alignedMinutiaDirection = CS.V[vm].alignedMinutiaDirection;

                                int indexTermInt = getIntFromBitArray(indexTermBits);

                                if (indexTable[idHashTable * Nd + k].ContainsKey(indexTermInt))
                                {
                                    indexTable[idHashTable * Nd + k][indexTermInt].Add(idxItem);
                                }
                                else
                                {
                                    List<IndexItem> listIdxItem = new List<IndexItem>();
                                    listIdxItem.Add(idxItem);
                                    indexTable[idHashTable * Nd + k].Add(indexTermInt, listIdxItem);
                                }
                            }
                        }
                        idHashTable++;
                    }
                }
            }
        }

        public void queryFPInIndexTableFromMCCText_AppendixB_WithPoseText(string queryMCCFolder, string queryPoseFolder, string saveResultFilePath,
                                                                          char enrollFIDCharPair = 'F', char queryFIDCharPair = 'S', 
                                                                          bool saveImgMatchedMinutiae = false, string saveImgMatchedMinutiaeFolder = "",
                                                                          string saveImgMatchedMinutiaeExtension = ".png", 
                                                                          string galleryMCCFolder = "", string galleryPoseFolder = "",
                                                                          string galleryFPImgFolder = "", string queryFPImgFolder = "",
                                                                          string gqImgExtension = ".png")
        {
            var vQueryMCCFileNames = new DirectoryInfo(queryMCCFolder).GetFileSystemInfos("*.txt").Select(x => x.FullName);
            string[] queryMCCFileNames = vQueryMCCFileNames.ToArray();

            string MCCFpath = Path.GetDirectoryName(queryMCCFileNames[0]);

            var vQueryPoseFileNames = new DirectoryInfo(queryPoseFolder).GetFileSystemInfos("*.txt").Select(x => x.FullName);
            string[] queryPoseFileNames = vQueryPoseFileNames.ToArray();

            string PoseFpath = Path.GetDirectoryName(queryPoseFileNames[0]);

            Stopwatch stopwatch = new Stopwatch();

            for (int f = 0; f < queryMCCFileNames.Length; f++)//107
            {
                // Begin timing.
                stopwatch.Restart();

                string fname = Path.GetFileNameWithoutExtension(queryMCCFileNames[f]);

                Pose pose = getPoseFromFile(PoseFpath + "//" + fname + ".txt");

                CylindersTemplate CS = new CylindersTemplate();
                CS.F_ID = fname;
                CS.pose = pose;
                createBitCylinderTemplateFromMCCFile_AppendixB(queryMCCFileNames[f], ref CS, pose);

                Dictionary<string, Dictionary<int, RetrieveItem>>[] mPairTable = new Dictionary<string, Dictionary<int, RetrieveItem>>[CS.V.Length];

                Dictionary<string, double> finalScoreTable = searchIndexTable(CS, mPairTable);
                
                string[] F_IDCandidatelist = finalScoreTable.Keys.ToArray();
                int rank = Array.IndexOf(F_IDCandidatelist, CS.F_ID.Replace(queryFIDCharPair, enrollFIDCharPair)) + 1;

                // Stop timing.
                stopwatch.Stop();

                allRankDb.Add(CS.F_ID, rank);
                saveRetrivalResults(saveResultFilePath, CS.F_ID, rank);

                if (saveImgMatchedMinutiae)
                {
                    Image<Bgr, Single> inputBgr = new Image<Bgr, float>(Path.Combine(queryFPImgFolder, fname + gqImgExtension));
                    if (rank != 0)
                    {
                        Image<Bgr, Single> mateFPBgr = new Image<Bgr, float>(Path.Combine(galleryFPImgFolder, F_IDCandidatelist[rank - 1] + gqImgExtension));

                        Pose gPose = getPoseFromFile(Path.Combine(galleryPoseFolder , F_IDCandidatelist[rank - 1] + ".txt"));
                        CylindersTemplate gCS = new CylindersTemplate();
                        gCS.F_ID = F_IDCandidatelist[rank - 1];
                        gCS.pose = gPose;
                        createBitCylinderTemplateFromMCCFile_AppendixB(Path.Combine(galleryMCCFolder, F_IDCandidatelist[rank - 1] + ".txt"), ref gCS, gPose);

                        Image<Bgr, Single> result = new Image<Bgr, float>(1, 1);
                        plotMatchedMinutiae(inputBgr, mateFPBgr, CS, gCS, mPairTable, ref result);
                        result.Save(Path.Combine(saveImgMatchedMinutiaeFolder, fname + "_r" + rank.ToString() + "_" + F_IDCandidatelist[rank - 1] + saveImgMatchedMinutiaeExtension));
                        result.Dispose();
                        mateFPBgr.Dispose();
                        
                        if (rank != 1)
                        {
                            Image<Bgr, Single> matchedRank1FPBgr = new Image<Bgr, float>(Path.Combine(galleryFPImgFolder, F_IDCandidatelist[0] + gqImgExtension));

                            gPose = getPoseFromFile(Path.Combine(galleryPoseFolder, F_IDCandidatelist[0] + ".txt"));
                            gCS = new CylindersTemplate();
                            gCS.F_ID = F_IDCandidatelist[0];
                            gCS.pose = gPose;
                            createBitCylinderTemplateFromMCCFile_AppendixB(Path.Combine(galleryMCCFolder, F_IDCandidatelist[0] + ".txt"), ref gCS, gPose);

                            Image<Bgr, Single> resultRank1 = new Image<Bgr, float>(1, 1);
                            plotMatchedMinutiae(inputBgr, matchedRank1FPBgr, CS, gCS, mPairTable, ref resultRank1);
                            resultRank1.Save(Path.Combine(saveImgMatchedMinutiaeFolder, fname + "_r1" + "_" + F_IDCandidatelist[0] + saveImgMatchedMinutiaeExtension));
                            resultRank1.Dispose();
                            matchedRank1FPBgr.Dispose();
                        }
                    }
                    inputBgr.Dispose();
                }

                finalScoreTable.Clear();
                finalScoreTable = null;
                for (int i = 0; i < mPairTable.Length; i++)
                {
                    if (mPairTable[i] != null)
                    {
                        var mPairTableKey = mPairTable[i].Keys.ToArray();
                        for (int j = 0; j < mPairTableKey.Length; j++)
                        {
                            mPairTable[i][mPairTableKey[j]].Clear();
                            mPairTable[i][mPairTableKey[j]] = null;
                        }
                        mPairTable[i].Clear();
                        mPairTable[i] = null;
                    }
                }
                mPairTable = null;
                Console.WriteLine("{0} is queried successfully. Time elapsed: {1}", fname, stopwatch.Elapsed);
            }

        }

        public Dictionary<string, double> searchIndexTable(CylindersTemplate CS, Dictionary<string, Dictionary<int, RetrieveItem>>[] mPairTable)
        {
            Dictionary<string, Dictionary<int, List<int>>> fIDMateTable = new Dictionary<string, Dictionary<int, List<int>>>();
            Dictionary<string, Dictionary<int, int>> countMIDForFID = new Dictionary<string, Dictionary<int, int>>();
            Dictionary<string, double> finalScoreTable = new Dictionary<string, double>();

            for (int mi = 0; mi < CS.V.Length; mi++)
            {
                if (CS.V[mi].valid)
                {
                    mPairTable[mi] = new Dictionary<string, Dictionary<int, RetrieveItem>>();

                    int idHashTable = 0;
                    foreach (int[][] hf in H)
                    {
                        for (int k = 0; k < Nd; k++)
                        {
                            BitArray indexTermBits = new BitArray(hashBits);

                            for (int hb = 0; hb < hashBits; hb++)
                            {
                                int getElementAt = hf[k][hb];
                                indexTermBits[hb] = CS.V[mi].cmBit[getElementAt];
                            }

                            if (popCount(indexTermBits) >= minPC)
                            {
                                int indexTermInt = getIntFromBitArray(indexTermBits);

                                if (indexTable[idHashTable * Nd + k].ContainsKey(indexTermInt))
                                {
                                    List<IndexItem> compatibleList = indexTable[idHashTable * Nd + k][indexTermInt].Where(
                                                                                        item =>
                                                                                        distanceOf2Points(item.MInFixedPoseCoordinate, CS.V[mi].MInFixedPoseCoordinate) <= eL
                                                                                        &&
                                                                                        Math.Abs(differenceAngleOf2Directions_0To2Pi(item.alignedMinutiaDirection, CS.V[mi].alignedMinutiaDirection)) <= eTheta
                                                                                        ).ToList();

                                    foreach (IndexItem item in compatibleList)
                                    {
                                        string F_ID = item.F_ID;
                                        int mj = item.M_ID;

                                        RetrieveItem rtItem = new RetrieveItem();
                                        rtItem.F_ID = item.F_ID;
                                        rtItem.M_ID = item.M_ID;
                                        rtItem.distanceFromPose = item.distanceFromPose;
                                        rtItem.angleFromPose = item.angleFromPose;
                                        rtItem.MInPoseCoordinate = new PointF(item.MInPoseCoordinate.X, item.MInPoseCoordinate.Y);
                                        rtItem.MInFixedPoseCoordinate = new PointF(item.MInFixedPoseCoordinate.X, item.MInFixedPoseCoordinate.Y);
                                        rtItem.alignedMinutiaDirection = item.alignedMinutiaDirection;
                                        rtItem.Hits = 1;

                                        if (mPairTable[mi].ContainsKey(F_ID))
                                        {
                                            if (mPairTable[mi][F_ID].ContainsKey(mj))
                                            {
                                                mPairTable[mi][F_ID][mj].Hits++;
                                            }
                                            else
                                            {
                                                int[] previousItem = mPairTable[mi][F_ID].Keys.ToArray();
                                                if (distanceOf2Points(mPairTable[mi][F_ID][previousItem[0]].MInFixedPoseCoordinate, CS.V[mi].MInFixedPoseCoordinate) >
                                                   distanceOf2Points(rtItem.MInFixedPoseCoordinate, CS.V[mi].MInFixedPoseCoordinate))
                                                {
                                                    mPairTable[mi][F_ID].Remove(previousItem[0]);
                                                    mPairTable[mi][F_ID].Add(mj, rtItem);
                                                }
                                            }
                                        }
                                        else
                                        {
                                            Dictionary<int, RetrieveItem> t = new Dictionary<int, RetrieveItem>();
                                            t.Add(mj, rtItem);
                                            mPairTable[mi].Add(F_ID, t);
                                        }
                                    }
                                    compatibleList.Clear();
                                    compatibleList.TrimExcess();
                                    compatibleList = null;
                                }
                                else
                                {
                                    //hashKey was not found in index table
                                }
                            }
                        }
                        idHashTable++;
                    }


                }
            }

            for (int mi = 0; mi < mPairTable.Length; mi++)
            {
                if (mPairTable[mi] != null)
                {
                    foreach (KeyValuePair<string, Dictionary<int, RetrieveItem>> fid_mj in mPairTable[mi])
                    {
                        string F_ID = fid_mj.Key;
                        List<int> listMj = fid_mj.Value.Keys.ToList();
                        foreach (int mj in listMj)
                        {
                            if (fIDMateTable.ContainsKey(F_ID))
                            {
                                if (fIDMateTable[F_ID].ContainsKey(mj))
                                {
                                    if (!fIDMateTable[F_ID][mj].Contains(mi))
                                    {
                                        fIDMateTable[F_ID][mj].Add(mi);
                                    }
                                }
                                else
                                {
                                    fIDMateTable[F_ID].Add(mj, new List<int>());
                                    fIDMateTable[F_ID][mj].Add(mi);
                                }
                            }
                            else
                            {
                                Dictionary<int, List<int>> t = new Dictionary<int, List<int>>();
                                t.Add(mj, new List<int>());
                                t[mj].Add(mi);
                                fIDMateTable.Add(F_ID, t);
                            }
                        }
                        listMj.Clear();
                        listMj.TrimExcess();
                        listMj = null;
                    }
                }
            }

            List<string> toRemoveFID = new List<string>();
            foreach (KeyValuePair<string, Dictionary<int, List<int>>> fid_mj in fIDMateTable)
            {
                List<int> toRemoveMj = new List<int>();
                foreach (KeyValuePair<int, List<int>> mjfid in fid_mj.Value)
                {
                    if (mjfid.Value.Count <= 1)
                    {
                        toRemoveMj.Add(mjfid.Key);
                    }
                }

                if (toRemoveMj.Count == fid_mj.Value.Count())
                {
                    toRemoveFID.Add(fid_mj.Key);
                }
                foreach (var key in toRemoveMj)
                {
                    fid_mj.Value.Remove(key);
                }

                toRemoveMj.Clear();
                toRemoveMj.TrimExcess();
                toRemoveMj = null;
            }

            foreach (var key in toRemoveFID)
            {
                fIDMateTable.Remove(key);
            }

            toRemoveFID.Clear();
            toRemoveFID.TrimExcess();
            toRemoveFID = null;

            foreach (KeyValuePair<string, Dictionary<int, List<int>>> fid_mj in fIDMateTable)
            {
                var duplicateMatchedMj = fid_mj.Value.SelectMany(c => c.Value).GroupBy(c => c).Where(c => c.Count() > 1).ToList();

                if (duplicateMatchedMj.Count > 0)
                {
                    //int error = 0;
                }
                else
                {
                    foreach (KeyValuePair<int, List<int>> mjfid in fid_mj.Value)
                    {
                        int Mj = mjfid.Key;
                        List<int> listMi = mjfid.Value;
                        double minDistance = 9999999;
                        int mateMi = -1;

                        for (int i = 0; i < listMi.Count; i++)
                        {
                            int Mi = listMi[i];
                            double distanceMiMj = distanceOf2Points(mPairTable[Mi][fid_mj.Key][Mj].MInFixedPoseCoordinate, CS.V[Mi].MInFixedPoseCoordinate);
                            if (distanceMiMj <= minDistance)
                            {
                                minDistance = distanceMiMj;
                                mateMi = Mi;
                            }
                        }

                        if (mateMi != -1)
                        {
                            List<int> toRemoveMi = listMi.Where(m => m != mateMi).ToList();

                            foreach (var keyMi in toRemoveMi)
                            {
                                mPairTable[keyMi].Remove(fid_mj.Key);
                            }

                            toRemoveMi.Clear();
                            toRemoveMi.TrimExcess();
                            toRemoveMi = null;
                        }
                        else
                        {
                            //int error = 1;
                        }

                        listMi.Clear();
                        listMi.TrimExcess();
                        listMi = null;
                    }
                }

                duplicateMatchedMj.Clear();
                duplicateMatchedMj.TrimExcess();
                duplicateMatchedMj = null;
            }

            //Sum Score
            for (int mi = 0; mi < mPairTable.Length; mi++)
            {
                if (mPairTable[mi] != null)
                {
                    foreach (KeyValuePair<string, Dictionary<int, RetrieveItem>> fid_mj in mPairTable[mi])
                    {
                        foreach (KeyValuePair<int, RetrieveItem> mjt in fid_mj.Value)
                        {
                            double t = mjt.Value.Hits;
                            double pp = p;
                            double hh = h;
                            double ph = pp / hh;
                            double Sigj = Math.Pow(t, ph);

                            if (finalScoreTable.ContainsKey(fid_mj.Key))
                            {
                                finalScoreTable[fid_mj.Key] += Sigj;
                            }
                            else
                            {
                                finalScoreTable.Add(fid_mj.Key, Sigj);
                            }
                        }
                    }
                }
            }

            foreach (string ng in finalScoreTable.Keys.ToArray())
            {
                finalScoreTable[ng] /= numOfMForFID[ng];
            }
            finalScoreTable = finalScoreTable.OrderByDescending(u => u.Value).ToDictionary(z => z.Key, y => y.Value);

            var fIDMateTableKey = fIDMateTable.Keys.ToArray();
            for (int ifid = 0; ifid < fIDMateTableKey.Length; ifid++)
            {
                string F_ID = fIDMateTableKey[ifid];
                var fIDMateTableKeyMjKey = fIDMateTable[F_ID].Keys.ToArray();

                for (int jfid = 0; jfid < fIDMateTableKeyMjKey.Length; jfid++)
                {
                    int mj = fIDMateTableKeyMjKey[jfid];
                    fIDMateTable[F_ID][mj].Clear();
                    fIDMateTable[F_ID][mj].TrimExcess();
                    fIDMateTable[F_ID][mj] = null;
                }
                fIDMateTable[F_ID].Clear();
                fIDMateTable[F_ID] = null;
            }
            fIDMateTable.Clear();
            fIDMateTable = null;
            
            foreach (KeyValuePair<string, Dictionary<int, int>> cid in countMIDForFID)
            {
                string F_ID = cid.Key;
                countMIDForFID[F_ID].Clear();
                countMIDForFID[F_ID] = null;
            }
            countMIDForFID.Clear();
            countMIDForFID = null;

            return finalScoreTable;
        }

        private int popCount(BitArray bitArray)
        {
            int count = 0;
            foreach (bool bit in bitArray)
            {
                if (bit)
                {
                    count++;
                }
            }
            return count;
        }

        private int getIntFromBitArray(BitArray bitArray)
        {
            if (bitArray.Length > 32)
                throw new ArgumentException("Argument length shall be at most 32 bits.");

            int[] array = new int[1];
            bitArray.CopyTo(array, 0);
            return array[0];
        }

        public void WriteIndexTableToBinary(string filePath)
        {
            using (FileStream fs = File.OpenWrite(filePath))
            using (BinaryWriter writer = new BinaryWriter(fs))
            {
                writer.Write(numOfMForFID.Count);
                foreach (var pair in numOfMForFID)
                {
                    writer.Write(pair.Key);
                    writer.Write(pair.Value);
                }

                // Put count.
                writer.Write(indexTable.Length);
                for (int i = 0; i < indexTable.Length; i++)
                {
                    writer.Write(indexTable[i].Count);
                    // Write pairs.
                    foreach (var pair in indexTable[i])
                    {
                        writer.Write(pair.Key);
                        writer.Write(pair.Value.Count);
                        for (int j = 0; j < pair.Value.Count; j++)
                        {
                            writer.Write(pair.Value[j].F_ID);
                            writer.Write(pair.Value[j].M_ID);
                            writer.Write(pair.Value[j].distanceFromPose);
                            writer.Write(pair.Value[j].angleFromPose);
                            writer.Write(pair.Value[j].MInPoseCoordinate.X);
                            writer.Write(pair.Value[j].MInPoseCoordinate.Y);
                            writer.Write(pair.Value[j].MInFixedPoseCoordinate.X);
                            writer.Write(pair.Value[j].MInFixedPoseCoordinate.Y);
                            writer.Write(pair.Value[j].alignedMinutiaDirection);
                        }
                    }
                }
            }
        }

        public void ReadIndexTableFromBinary(string filePath)
        {
            var result = new List<Dictionary<int, List<IndexItem>>>();
            using (FileStream fs = File.OpenRead(filePath))
            using (BinaryReader reader = new BinaryReader(fs))
            {
                // Get count.
                int numOfMForFIDCount = reader.ReadInt32();
                for (int i = 0; i < numOfMForFIDCount; i++)
                {
                    string F_ID = reader.ReadString();
                    int numOfM = reader.ReadInt32();
                    numOfMForFID.Add(F_ID, numOfM);
                }

                // Get count.
                int count = reader.ReadInt32();
                // Read in all pairs.
                for (int i = 0; i < count; i++)
                {
                    Dictionary<int, List<IndexItem>> dicTmp = new Dictionary<int, List<IndexItem>>();
                    int numOfKey = reader.ReadInt32();
                    for (int j = 0; j < numOfKey; j++)
                    {
                        int key = reader.ReadInt32();
                        int numOfItem = reader.ReadInt32();
                        List<IndexItem> listIdxItem = new List<IndexItem>();
                        for (int k = 0; k < numOfItem; k++)
                        {
                            IndexItem idxItem = new IndexItem();
                            idxItem.F_ID = reader.ReadString();
                            idxItem.M_ID = reader.ReadInt32();
                            idxItem.distanceFromPose = reader.ReadDouble();
                            idxItem.angleFromPose = reader.ReadDouble();
                            float x = reader.ReadSingle();
                            float y = reader.ReadSingle();
                            idxItem.MInPoseCoordinate = new PointF(x, y);
                            x = reader.ReadSingle();
                            y = reader.ReadSingle();
                            idxItem.MInFixedPoseCoordinate = new PointF(x, y);
                            idxItem.alignedMinutiaDirection = reader.ReadDouble();
                            listIdxItem.Add(idxItem);
                        }

                        dicTmp.Add(key, listIdxItem);
                    }

                    result.Add(dicTmp);
                }
            }
            indexTable = result.ToArray();
        }

        private void saveRetrivalResults(string filePath, string F_ID, int rank)
        {
            using (var writer = new StreamWriter(filePath, true))
            {
                writer.WriteLine("{0}\t{1}", F_ID, rank);
            }
        }
        #endregion

        #region Utility
        public int calLin(int i, int j, int k, int Ns)
        {
            return ((k - 1) * (Ns * Ns)) + ((j - 1) * Ns) + i;
        }

        public double distanceOf2Points(PointF a, PointF b)
        {
            return Math.Sqrt(Math.Pow((a.X - b.X), 2) + Math.Pow(a.Y - b.Y, 2));
        }

        public double angleOf2Points(PointF a, PointF b)
        {
            return Math.Atan2(a.Y - b.Y, a.X - b.X);
        }

        /// <summary>
        /// Calculate diffence between two angles. The result is in range [-Pi, Pi]
        /// </summary>
        /// <param name="theta1">input angle1 in range 0 to 2Pi</param>
        /// <param name="theta2">input angle2 in range 0 to 2Pi</param>
        /// <returns>
        /// Return diffence of theta1 - theta2 if in range [-Pi, Pi] 
        /// or return 2Pi + diffence if the diffence less than -Pi,
        /// or return -2Pi + diffence if the diffence greater than Pi.
        /// </returns>
        public double differenceAngleOf2Directions_0To2Pi(double theta1, double theta2)
        {
            double result = theta1 - theta2;
            if (result < -Math.PI)
            {
                result = (2.0 * Math.PI) + result;
            }
            else if (result >= Math.PI)
            {
                result = (-2.0 * Math.PI) + result;
            }
            return result;
        }
        #endregion

        #region Plot
        public void plotPoseOnImage(Pose p, ref Image<Bgr, Single> inputBgr)
        {
            PointF str1 = new PointF(p.p.X, p.p.Y);
            PointF end1 = new PointF();
            PointF endArrow1 = new PointF();
            PointF endArrow2 = new PointF();
            float length = 200;
            float lengthArrow = 30;
            end1.X = (float)(str1.X + length * Math.Cos(p.angle));
            end1.Y = (float)(str1.Y + length * Math.Sin(p.angle));
            endArrow1.X = (float)(end1.X - lengthArrow * Math.Cos(p.angle + Math.PI / 7));
            endArrow1.Y = (float)(end1.Y - lengthArrow * Math.Sin(p.angle + Math.PI / 7));
            endArrow2.X = (float)(end1.X - lengthArrow * Math.Cos(p.angle - Math.PI / 7));
            endArrow2.Y = (float)(end1.Y - lengthArrow * Math.Sin(p.angle - Math.PI / 7));

            LineSegment2DF poseLine = new LineSegment2DF(str1, end1);
            LineSegment2DF arrowLine1 = new LineSegment2DF(end1, endArrow1);
            LineSegment2DF arrowLine2 = new LineSegment2DF(end1, endArrow2);

            inputBgr.Draw(new CircleF(p.p, 7), new Bgr(0, 0, 255), 3);
            inputBgr.Draw(poseLine, new Bgr(0, 0, 255), 3);
            inputBgr.Draw(arrowLine1, new Bgr(0, 0, 255), 3);
            inputBgr.Draw(arrowLine2, new Bgr(0, 0, 255), 3);

            //ImageViewer.Show(inputBgr);
        }

        public void plotMinutiaeFromCylinderTemplate(Image<Bgr, Single> FPImage, CylindersTemplate template)
        {
            for (int i = 0; i < template.V.Length; i++)
            {
                double alignedMinutiaDirection = template.V[i].alignedMinutiaDirection;
                float alignedMX = template.V[i].MInPoseCoordinate.X;
                float alignedMY = template.V[i].MInPoseCoordinate.Y;
                float length = 15;
                PointF str = new PointF(alignedMX, alignedMY);
                PointF end = new PointF();
                end.X = (float)(str.X + length * Math.Cos(alignedMinutiaDirection));
                end.Y = (float)(str.Y + length * Math.Sin(alignedMinutiaDirection));
                LineSegment2DF alignedMinutiaLine = new LineSegment2DF(str, end);
                FPImage.Draw(new CircleF(str, 5), new Bgr(0, 0, 255));
                FPImage.Draw(alignedMinutiaLine, new Bgr(0, 0, 255), 1);

                Point mIDPoint = new Point();
                mIDPoint.X = (int)str.X - 8;
                mIDPoint.Y = alignedMinutiaDirection > Math.PI ? (int)str.Y + 13 : (int)str.Y - 8;
                FPImage.Draw(template.V[i].M_ID.ToString(), mIDPoint, FontFace.HersheyPlain, 0.75, new Bgr(Color.Red));
            }
        }

        public void plotMatchedMinutiae(Image<Bgr, Single> queryFPImage, Image<Bgr, Single> galleryFPImage, CylindersTemplate queryTemplate, CylindersTemplate galleryTemplate, Dictionary<string, Dictionary<int, RetrieveItem>>[] mPairTable, ref Image<Bgr, Single> result)
        {
            Image<Bgr, Single> imgQueryAlignedMinutiaeWithPose = queryFPImage.Rotate(queryTemplate.rotateDirectionDegree, new Bgr(Color.Black)).Clone();
            plotPoseOnImage(queryTemplate.alignedPose, ref imgQueryAlignedMinutiaeWithPose);

            Image<Bgr, Single> imgGalleryAlignedMinutiaeWithPose = galleryFPImage.Rotate(galleryTemplate.rotateDirectionDegree, new Bgr(Color.Black)).Clone();
            plotPoseOnImage(galleryTemplate.alignedPose, ref imgGalleryAlignedMinutiaeWithPose);

            plotMinutiaeFromCylinderTemplate(imgQueryAlignedMinutiaeWithPose, queryTemplate);
            plotMinutiaeFromCylinderTemplate(imgGalleryAlignedMinutiaeWithPose, galleryTemplate);

            Image<Bgr, Single> imgMatched = imgQueryAlignedMinutiaeWithPose.ConcateHorizontal(imgGalleryAlignedMinutiaeWithPose).Clone();
            Image<Bgr, Single> imgFixedPoseMatched = new Image<Bgr, float>(queryFPImage.Width, queryFPImage.Height, new Bgr(Color.White));
            Pose fixedPose = new Pose();
            fixedPose.p = new PointF((float)(queryFPImage.Width / 2.0), (float)(queryFPImage.Height / 2.0));
            fixedPose.angle = -(Math.PI / 2.0);
            plotPoseOnImage(fixedPose, ref imgFixedPoseMatched);

            for (int mi = 0; mi < mPairTable.Length; mi++)
            {
                if (mPairTable[mi] != null)
                {
                    var mate = mPairTable[mi].Where(i => i.Key == galleryTemplate.F_ID);
                    foreach (KeyValuePair<string, Dictionary<int, RetrieveItem>> fid_mj in mate)
                    {
                        foreach (KeyValuePair<int, RetrieveItem> mjt in fid_mj.Value)
                        {
                            PointF miP = queryTemplate.V[mi].MInPoseCoordinate;
                            PointF mjP = mjt.Value.MInPoseCoordinate;
                            mjP.X += imgQueryAlignedMinutiaeWithPose.Width;
                            LineSegment2DF alignedMinutiaLine = new LineSegment2DF(miP, mjP);
                            imgMatched.Draw(mjt.Value.Hits.ToString(), Point.Round(miP), FontFace.HersheyComplexSmall, 1.0, new Bgr(Color.Blue));
                            imgMatched.Draw(alignedMinutiaLine, new Bgr(255, 0, 0), 1);
                            //ImageViewer.Show(imgMatched);

                            float mLength = 15;
                            PointF miPStr = queryTemplate.V[mi].MInFixedPoseCoordinate;
                            PointF miPEnd = new PointF();
                            miPEnd.X = (float)(miPStr.X + mLength * Math.Cos(queryTemplate.V[mi].alignedMinutiaDirection));
                            miPEnd.Y = (float)(miPStr.Y + mLength * Math.Sin(queryTemplate.V[mi].alignedMinutiaDirection));
                            LineSegment2DF miLine = new LineSegment2DF(miPStr, miPEnd);
                            imgFixedPoseMatched.Draw(new CircleF(miPStr, 5), new Bgr(0, 0, 255));
                            imgFixedPoseMatched.Draw(miLine, new Bgr(0, 0, 255), 1);

                            PointF mjPStr = mjt.Value.MInFixedPoseCoordinate;
                            PointF mjPEnd = new PointF();
                            mjPEnd.X = (float)(mjPStr.X + mLength * Math.Cos(mjt.Value.alignedMinutiaDirection));
                            mjPEnd.Y = (float)(mjPStr.Y + mLength * Math.Sin(mjt.Value.alignedMinutiaDirection));
                            LineSegment2DF mjLine = new LineSegment2DF(mjPStr, mjPEnd);
                            imgFixedPoseMatched.Draw(new CircleF(mjPStr, 5), new Bgr(0, 255, 0));
                            imgFixedPoseMatched.Draw(mjLine, new Bgr(0, 255, 0), 1);

                            LineSegment2DF connectedMiMjLine = new LineSegment2DF(miPStr, mjPStr);
                            imgFixedPoseMatched.Draw(connectedMiMjLine, new Bgr(255, 0, 0), 1);
                            //ImageViewer.Show(imgMatched.ConcateHorizontal(imgFixedPoseMatched));
                        }
                    }
                }
            }
            //ImageViewer.Show(imgMatched.ConcateHorizontal(imgFixedPoseMatched));

            result = imgMatched.ConcateHorizontal(imgFixedPoseMatched).Clone();

            imgQueryAlignedMinutiaeWithPose.Dispose();
            imgGalleryAlignedMinutiaeWithPose.Dispose();
            imgMatched.Dispose();
            imgFixedPoseMatched.Dispose();
        }
        #endregion
    }
}
