# Fingerprint Indexing with Pose Constraint (Yijing Su, Jianjiang Feng, and Jie Zhou, 2016)

This repository provides an implementation of the fingerprint indexing with pose constraints in C# using EmguCV. The concept was proposed by Yijing Su, Jianjiang Feng, and Jie Zhou, 2016. Please refer to [Su2016](https://www.sciencedirect.com/science/article/abs/pii/S003132031600008X) for more details.

### Please note that the source code is unofficial and is implemented for study purpose only. The correctness has not been confirmed by the author.

The source code is organized as one class for integration flexibility.

#### Usage
Open the solution with Microsoft Visual Studio and please 'Rebuild' the solution to gather all required packages.

#### Contents
1. [Requirement](#requirement)
2. [Input and Output Format](#input-and-output-format)
3. [Demo](#demo)

#### Requirement
EmguCV 4.3, .NET Framework 4.0

#### Input and Output Format

##### Input Format
This approach requires two features as input as follows:
	
  1. ***Mcc Template File Text Format***, please refer to [MccSdk Documentation v2.0](http://biolab.csr.unibo.it/research.asp?organize=Activities&select=&selObj=82&pathSubj=111%7C%7C8%7C%7C82&Req=&), Appendix B. This work converts minutiae format from conventional algorithm to **Minutiae Template File Text Format**, Appendix A, then converts to ***Mcc Template File Text Format*** using MccSdk V2.0 with MccIndexPaperEnrollParameters.xml.
  
  2. Pose text file contains position (x, y) and direction (θ), where θ=[-90°,90°], 0° points in upward direction and increases in counter-clockwise direction. The formatting syntax is ```x y θ```

##### Output Format

The output has been provided as a text file which contains a query fingerprint ID (filename) and follows by retrieved rank, i.e., ```Filename Rank```. For example, ```S0024301 1```. Note that the result file will be written in append mode.

Furthermore, matched minutiae of query and gallery fingerprints can be saved as image file. 
If the pair fingerprint is retrieved at rank-1, there will be one matched minutiae image. 
Otherwise, there will be two matched minutiae images. The first image is for the rank-1 retrieved gallery fingerprint and another is for the rank-k retrieved pair fingerprint. 

The image filename format is ```QueryFilename_rRank_RetrivedFilename.xxx```. For example, ``` S0024301_r1_F0024301.png```. The number displayed next to each minutia on the left image indicates the number of hash functions that the minutia on the lefe image collides with a minutia on the middle image.

![Example result](/Dataset/Result/SExample001_r1_FExample001.png)

#### Demo

There are 3 steps to run fingerprint indexing with pose constraints.

1. Create hash functions:
	
	The hash functions can be created by calling construction function of ```FPIndexingMCCPoseConst_Su2016``` class as follows.
	
	* To create new hash functions
	```c#
	FPIndexingMCCPoseConst_Su2016 su2016 = new FPIndexingMCCPoseConst_Su2016();
	```
	* To load hash functions from text file
	```c#
	string projectDirectory = Directory.GetParent(Environment.CurrentDirectory).Parent.Parent.FullName;
	
	FPIndexingMCCPoseConst_Su2016 su2016 = new FPIndexingMCCPoseConst_Su2016(projectDirectory + @"\Dataset\Hash.txt");
	```
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;or
	```c#
	FPIndexingMCCPoseConst_Su2016 su2016 = new FPIndexingMCCPoseConst_Su2016(projectDirectory + @"\Dataset\Hash.txt", 60.0, (10.0 * Math.PI) / 180.0);
	```
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;to specify location uncertainty (eL) and direction uncertainty (eTheta) parameters (please refer to [Su2016](https://www.sciencedirect.com/science/article/abs/pii/S003132031600008X) for more details). This project sets eL = 80.0 and eTheta = (20.0 * Math.PI) / 180.0 as default value. If readers use pose from [Ouyang2017](https://ieeexplore.ieee.org/document/8272707) eL and eTheta should be adjusted to 60.0 and (10.0 * Math.PI) / 180.0, respectively.

    After the hash functions are created, they can be saved to text file for further usage by calling ```writeHashFunctionsToText``` function as follows.
    ```c#
    su2016.writeHashFunctionsToText(projectDirectory + @"\Dataset\Hash.txt");
    ```
    
2. Enroll gallery fingerprints to index table:
	
	This step transforms minutiae position and direction into pose coordinate system. 
	
	Then it generates the index terms for each minutia by using MCC features and hash functions. 
	
	Finally, it enrolls those index terms with fingerprint ID (filename) and transformed minutiae position and direction.
	
	This step can be done by calling ```enrollFPFromMCCText_AppendixB_WithPoseText``` function as follows.
	```c#
	su2016.enrollFPFromMCCText_AppendixB_WithPoseText(projectDirectory + @"\Dataset\MCCTemplateAppendixB\Gallery",
                                                          projectDirectory + @"\Dataset\Pose");
	```
    
    Where ```projectDirectory + @"\Dataset\MCCTemplateAppendixB\Gallery"``` specifies the gallery ***Mcc Template File Text Format*** path and ```projectDirectory + @"\Dataset\Pose"``` specifies the pose text file path.
    
    After the index table is created, it can be saved to binary file for further usage by calling ```WriteIndexTableToBinary``` function as follows.
    ```c#
    su2016.WriteIndexTableToBinary(projectDirectory + @"\Dataset\IndexTable.bin");
    ```
    And the ```ReadIndexTableFromBinary``` function can be called to load the index table binary file as follows.
    ```c#
    su2016.ReadIndexTableFromBinary(projectDirectory + @"\Dataset\IndexTable.bin");
    ```
	
3. Search query fingerprints from index table:
	
	This step retrieves gallery fingerprints from index table by using index terms of each query fingerprint. The search of query fingerprints is performed by calling ```queryFPInIndexTableFromMCCText_AppendixB_WithPoseText``` function as follows.
	
	```c#
	su2016.queryFPInIndexTableFromMCCText_AppendixB_WithPoseText(projectDirectory + @"\Dataset\MCCTemplateAppendixB\Query",
                                                                     projectDirectory + @"\Dataset\Pose",
                                                                     projectDirectory + @"\Dataset\Result\result.txt",
                                                                     'F', 'S', true,
                                                                     projectDirectory + @"\Dataset\Result", ".png",
                                                                     projectDirectory + @"\Dataset\MCCTemplateAppendixB\Gallery",
                                                                     projectDirectory + @"\Dataset\Pose",
                                                                     projectDirectory + @"\Dataset\FPImage", projectDirectory + @"\Dataset\FPImage", ".bmp");
	```
    	
								     

    The details of each parameter are described as follows.
    1. ```projectDirectory + @"\Dataset\MCCTemplateAppendixB\Query"``` specifies the query ***Mcc Template File Text Format*** path
    2. ```projectDirectory + @"\Dataset\Pose"``` specifies the pose text file path of the query fingerprints
    3. ```projectDirectory + @"\Dataset\Result\result.txt"``` specifies the result text file of a query fingerprint ID (filename) and follows by retrieved rank
    4. ```'F'``` specifies the letter of gallery fingerprint ID (filename) in index table. By replacing this letter with the specific letter of query fingerprint, the retrieved fingerprint ID will be the same as query fingerprint ID and the pair gallery fingerprint can be found. The default value is set to ```'F'```.
    5. ```'S'``` specifies the letter of query fingerprint ID (filename). The default value is set to ```'S'```.
    6. ```true``` this parameter sets to true for saving retrieval result images as aforementioned in **Output Format** section. The default value is set to ```false```.
    7. ```projectDirectory + @"\Dataset\Result"``` specifies the path that the retrieval result images will be saved. The default value is set to ```""```.
    8. ```".png"``` specifies the extension of the retrieval result images. The default value is set to ```".png"```.
    9. ```projectDirectory + @"\Dataset\MCCTemplateAppendixB\Gallery"``` specifies the gallery ***Mcc Template File Text Format*** path for drawing the retrieval result images. The default value is set to ```""```.
    10. ```projectDirectory + @"\Dataset\Pose"``` specifies the pose text file path of the gallery fingerprints for drawing the retrieval result images. The default value is set to ```""```.
    11. ```projectDirectory + @"\Dataset\FPImage"``` specifies the gallery fingerprint images file path for drawing the retrieval result images. The default value is set to ```""```.
    12. ```projectDirectory + @"\Dataset\FPImage"``` specifies the query fingerprint images file path for drawing the retrieval result images. The default value is set to ```""```.
    13. ```".png"``` specifies the extension of the gallery and query fingerprint images. The default value is set to ```".png"```.


The complete step should be similar as follows.
```c#
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
```









