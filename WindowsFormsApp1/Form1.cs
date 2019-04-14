using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using System.Drawing.Imaging;
using System.Text.RegularExpressions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace WindowsFormsApp1
{
   
    public partial class CR : Form
    {
        double[,] forest1 = { { 0.7, 0.9, 0.36, 0.4 },{0.7,1.1,0.4,0.8 },{ 0.7,0.9,0.35,0.5} };
        double[] city2 = { -25.0088, -23, 3300, 839.4985, 933.3491 };
        double[] river3 = { -2.2365, -1.0706, 29.2263, 65.0149 };
        public CR()
        {
            InitializeComponent();
        }
        //Рекурсия
        public void Identification(Color[,] matr, int k)
        {
            int k0, k1, k2, k3;
            int k0c, k1c, k2c, k3c;
            int p = matr.GetLength(1);
            int n = 1;
            while ((k - System.Math.Pow(4, n)) >= 0)
            {
                k = k - (int)System.Math.Pow(4, n);
                n++;
            }

            int s = (int)System.Math.Pow(2, n);
            k0 = 4 * s * (k / s) + 2 * (k % s);
            k1 = 4 * s * (k / s) + 2 * (k % s) + 1;
            k2 = 4 * s * (k / s) + 2 * (k % s) + 2 * s;
            k3 = 4 * s * (k / s) + 2 * (k % s) + 2 * s + 1;
            if (n < 5)
            {
                while (n >= 1)
            {
                k0 += (int)System.Math.Pow(4, n);
                k1 += (int)System.Math.Pow(4, n);
                k2 += (int)System.Math.Pow(4, n);
                k3 += (int)System.Math.Pow(4, n);
                n--;
            }

            k0c = Eq(Kf(matr, k0));
            k1c = Eq(Kf(matr, k1));
            k2c = Eq(Kf(matr, k2));
            k3c = Eq(Kf(matr, k3));
            
                if (k0c != 0)
                    Pnt(matr, k0, k0c);
                else Identification(matr, k0);
                if (k1c != 0)
                    Pnt(matr, k1, k1c);
                else Identification(matr, k1);
                if (k2c != 0)
                    Pnt(matr, k2, k2c);
                else Identification(matr, k2);
                if (k3c != 0)
                    Pnt(matr, k3, k3c);
                else Identification(matr, k3);
            }
            else
            {
                while (n >=1)
                {
                    k0 += (int)System.Math.Pow(4, n);
                    k1 += (int)System.Math.Pow(4, n);
                    k2 += (int)System.Math.Pow(4, n);
                    k3 += (int)System.Math.Pow(4, n);
                    n--;
                }

                k0c = Eq(Kf(matr, k0));
                k1c = Eq(Kf(matr, k1));
                k2c = Eq(Kf(matr, k2));
                k3c = Eq(Kf(matr, k3));
                Pnt(matr, k0, k0c);
                Pnt(matr, k1, k1c);
                Pnt(matr, k2, k2c);
                Pnt(matr, k3, k3c);
                return;
            }
        }
        //Начальная
        public void Identification(Color[,] matr)
        {
            int k0, k1, k2, k3;
            int k0c, k1c, k2c, k3c;
            k0 = 0; k1 = 1; k2 = 2; k3 = 3;
            textBox1.Text = 1.ToString();
            k0c = Eq(Kf(matr, k0));
            k1c = Eq(Kf(matr, k1));
            k2c = Eq(Kf(matr, k2));
            k3c = Eq(Kf(matr, k3));
            if (k0c != 0)
                Pnt(matr, k0, k0c);
            else Identification(matr, k0);
            if (k1c != 0)
                Pnt(matr, k1, k1c);
            else Identification(matr, k1);
            if (k2c != 0)
                Pnt(matr, k2, k2c);
            else Identification(matr, k2);
            if (k3c != 0)
                Pnt(matr, k3, k3c);
            else Identification(matr, k3);

            return;
        }
        //поиск коэффициентов образа
        public double[,] Kf(Color [,] pix, int k)
        {
            
            int p = pix.GetLength(0);
            int n = 1;
            

            while ((k - System.Math.Pow(4, n))>=0)
            {
                k = k - (int)System.Math.Pow(4, n);
                n++;
            }
            
            textBox1.Text = n.ToString();
            int s = (int)System.Math.Pow(2, n);
            Color[,] m = new Color[p/s, p/s];
            for (int i = (k / s) * p / s , ii=0; i < ((k / s) + 1) * p / s; i++, ii++)
                for (int j = (k % s) * p / s , jj=0; j < ((k % s) + 1) * p / s; j++ , jj++)
                    m[ii,jj]= pix[i, j] ;
            double[,] arr = new double[3,2];
            ///////////////////////
            double[,] R = new double[p / s, p / s];
            double[,] G = new double[p / s, p / s];
            double[,] B = new double[p / s, p / s];
            for (int y = 0; y < p/s; y++)
                for (int x = 0; x < p/s; x++)
                {
                    R[x, y] = m[x, y].R;
                    G[x, y] = m[x, y].G;
                    B[x, y] = m[x, y].B;

                }

           //R//
            var mtr = Matrix<double>.Build.DenseOfArray(R);
            var svd = mtr.Svd(true);//разложение изображения
            double[] sing = new double[pix.GetLength(0)];
            for (int i = 1 ; i <  p / s; i++)
                sing[i] = svd.W[i, i];
            double sumx = 0, sumx1 = 0, sumy = 0, sumxy = 0;
            for (int i = 1; i <p/s; i++)
            {
                sumx += i;
                sumy += sing[i];
                sumx1 += i * i;
                sumxy += i * sing[i];
            }
            arr[0,0] = ((sing.Length) * sumxy - (sumx * sumy)) / ((sing.Length) * sumx1 - sumx * sumx);
            arr[0,1] = (sumy - arr[0,0] * sumx) / ((sing.Length));
           //G//
             mtr = Matrix<double>.Build.DenseOfArray(G);
             svd = mtr.Svd(true);//разложение изображения
             sing = new double[pix.GetLength(0)];
            for (int i = 1; i < p / s; i++)
                sing[i] = svd.W[i, i];
            sumx = 0; sumx1 = 0; sumy = 0; sumxy = 0;
            for (int i = 1; i < p / s; i++)
            {
                sumx += i;
                sumy += sing[i];
                sumx1 += i * i;
                sumxy += i * sing[i];
            }
            arr[1, 0] = ((sing.Length) * sumxy - (sumx * sumy)) / ((sing.Length) * sumx1 - sumx * sumx);
            arr[1, 1] = (sumy - arr[1, 0] * sumx) / ((sing.Length));

           //B//
            mtr = Matrix<double>.Build.DenseOfArray(B);
            svd = mtr.Svd(true);//разложение изображения
            sing = new double[pix.GetLength(0)];
            for (int i = 1; i < p / s; i++)
                sing[i] = svd.W[i, i];
            sumx = 0; sumx1 = 0; sumy = 0; sumxy = 0;
            for (int i = 1; i < p / s; i++)
            {
                sumx += i;
                sumy += sing[i];
                sumx1 += i * i;
                sumxy += i * sing[i];
            }
            arr[2, 0] = ((sing.Length) * sumxy - (sumx * sumy)) / ((sing.Length) * sumx1 - sumx * sumx);
            arr[2, 1] = (sumy - arr[2, 0] * sumx) / ((sing.Length));



            if ((n == 5) &&
                 (
                  (k >= 0 && k < 95) || (k >= 128 && k < 133) ||
                  (k / s >= 4 && k / s <= 8 && k % s >= 26 && k % s <= 30)||
                  (k / s >= 29 && k / s <= 31 && k % s >= 20 && k % s <= 23)|| (k / s >= 24 && k / s <= 26 && k % s >= 17 && k % s <= 18)|| k==919||k==918||k==755||k==787||k==819
                ) 
              )
            {
                using (StreamWriter svdfile = new StreamWriter(@"futest.txt", append: true))
                {
                    // svdfile.Write(i);
                    svdfile.Write(k); svdfile.Write(",");                 
                    svdfile.Write(arr[0,0]); svdfile.Write(",");
                    svdfile.Write(arr[0,1]); svdfile.Write(",");

                    svdfile.Write(arr[1, 0]); svdfile.Write(",");
                    svdfile.Write(arr[1, 1]); svdfile.Write(",");

                    svdfile.Write(arr[2, 0]); svdfile.Write(",");
                    svdfile.Write(arr[2, 1]); svdfile.WriteLine();
                    Color newColor = Color.FromArgb(0, 0, 0);
                    int i1 = (k / s) * p / s;
                    int i2 = ((k / s) + 1) * p / s-1;
                    int j1 = (k % s) * p / s;
                    int j2 = ((k % s) + 1) * p / s-1;

                    
                     for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                     {

                         pix[i1, j] = newColor;
                         pix[i1 + 1, j] = newColor;
                     }

                     for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                     {
                         pix[i, j1] = newColor;
                         pix[i, j2] = newColor;
                         pix[i, j1 + 1] = newColor;
                         pix[i, j2 - 1] = newColor;
                     }
                     for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                     {
                         pix[i2, j] = newColor;
                         pix[i2 -1, j] = newColor;
                     }
                };
            }
            
            return arr;
        }
        //Поиск для полного
        public double[] Kf(double[,] pix)
        {
            
            double[] arr = new double[2];
            var mtr = Matrix<double>.Build.DenseOfArray(pix);
            var svd = mtr.Svd(true);//разложение изображения
            double[] sing = new double[pix.GetLength(0)];
            for (int i = 1; i <16; i++)
                sing[i] = svd.W[i, i];
            double sumx = 0, sumx1 = 0, sumy = 0, sumxy = 0;
            for (int i = 1; i <16; i++)
            {
                sumx += i;
                sumy += sing[i];
                sumx1 += i * i;
                sumxy += i * sing[i];
            }
            arr[0] = ((sing.Length) * sumxy - (sumx * sumy)) / ((sing.Length) * sumx1 - sumx * sumx);
            arr[1] = (sumy - arr[0] * sumx) / ((sing.Length));
            return arr;

        }

        //Соответствие
        public int Eq(double[,] kf)
        {
            int n=0;
            /* if (kf[0] == 0 && kf[1] == 0)
                 n = 1;*/
            if ((kf[0, 0] >= forest1[0, 0]) && (kf[0, 0] <= forest1[0, 1]) && (kf[0, 1] >= forest1[0, 2]) && (kf[0, 1] <= forest1[0, 3])&&
                (kf[1, 0] >= forest1[1, 0]) && (kf[1, 0] <= forest1[1, 1]) && (kf[1, 1] >= forest1[1, 2]) && (kf[1, 1] <= forest1[1, 3])&&
                (kf[2, 0] >= forest1[2, 0]) && (kf[2, 0] <= forest1[2, 1]) && (kf[2, 1] >= forest1[2, 2]) && (kf[2, 1] <= forest1[2, 3]))
                n = 1;
            else if ((kf[1,0] >= city2[0]) && (kf[1,0] <= city2[1]) && (kf[1,1] >= city2[2]) && (kf[1,1] <= city2[3]))
                n=2;
            else if ((kf[2,0] >= river3[0]) && (kf[2,0] <= river3[1]) && (kf[2,1] >= river3[2]) && (kf[2,1] <= river3[3]))
                n=3;
            return n;

        }
        //Закрашивание
        public void Pnt(Color [,] pix, int k,int c)
        {
            int p = pix.GetLength(1);
            int n = 1;
            Color newColor=Color.FromArgb(0,0,0) ;
            while ((k - System.Math.Pow(4, n)) >= 0)
            {
                k = k - (int)System.Math.Pow(4, n);
                n++;
            }
            int s = (int)System.Math.Pow(2, n);


            /*switch (c)
            {
                case 1:
                    newColor = Color.FromArgb(0, 250, 0);
                    for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                        for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                            pix[i, j] = newColor;
                    break;
                case 2:
                    newColor = Color.FromArgb(250, 0, 0);
                    for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                        for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                            pix[i, j] = newColor;
                    break;
                case 3:
                    newColor = Color.FromArgb(0, 0, 250);
                    for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                        for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                            pix[i, j] = newColor;
                    break;
                default:
                   /* newColor = Color.FromArgb(0, 0, 0);
                    for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                        for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                            pix[i, j] = newColor;
                    break;
            }*/
            


        }
        //Закрашивание через файл
        public void Pnt16(Color[,] pix, int k, int c)
        {
            int p = pix.GetLength(1);
            int n = 5;
            Color newColor = Color.FromArgb(0, 0, 0);
           
            int s = (int)System.Math.Pow(2, n);


            switch (c)
            {
                case 0:
                    newColor = Color.FromArgb(0, 250, 0);
                    for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                        for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                            pix[i, j] = newColor;
                    break;
                case 1:
                    newColor = Color.FromArgb(250, 0, 0);
                    for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                        for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                            pix[i, j] = newColor;
                    break;
                case 2:
                    newColor = Color.FromArgb(0, 0, 250);
                    for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                        for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                            pix[i, j] = newColor;
                    break;
                default:
                    /*newColor = Color.FromArgb(0, 0, 0);
                    for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                        for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                            pix[i, j] = newColor;*/
                    break;
            }



        }
        //Основная часть
        private void button1_Click(object sender, EventArgs e)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
            Bitmap bmp = (Bitmap)pictureBox2.Image;//создаём объект битмапа для записи в него чёрнобелого изображения
            int M = bmp.Height;
            int N = bmp.Width;
            Color[,] color = new Color[bmp.Width, bmp.Height];
            for (int x = 0; x < bmp.Height; x++)
                for (int y = 0; y < bmp.Width; y++)
                {
                    color[y, x] = bmp.GetPixel(x, y);
                }







            Identification(color);


            string file = File.ReadAllText("prdct_knn.txt");
            int[] buf = file
                .Split(new char[] { ' '}, StringSplitOptions.RemoveEmptyEntries)
                .Select(n => int.Parse(n))
                .ToArray();



           /* for (int i = 0; i < 1024; i++)
                Pnt16(color, i, buf[i]);*/

                
               
           
            textBox1.Text = "end id";
            Bitmap image1 = new Bitmap(bmp.Height, bmp.Width);
            // Loop through the images pixels to reset color.
            for (int x = 0; x < image1.Height; x++)
            {
                for (int y = 0; y < image1.Width; y++)
                { 
                    image1.SetPixel(x, y, color[y,x]);
                }
            }

            // Set the PictureBox to display the image.
            pictureBox3.Image = image1;
            return;


        }

        //сохранение результата работы
        private void saveToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Bitmap bmpSave = (Bitmap)pictureBox3.Image;
            SaveFileDialog sfd = new SaveFileDialog();
            sfd.DefaultExt = "bmp";
            sfd.Filter = "Image files (*.bmp)|*.bmp|All files (*.*)|*.*";
            if (sfd.ShowDialog() == DialogResult.OK)
                bmpSave.Save(sfd.FileName, ImageFormat.Bmp);

        }
        //обесцвечивание входнового изображения
        private void openToolStripMenuItem_Click(object sender, EventArgs e)
        {
            // Диалог открытия файла.
            OpenFileDialog ofd = new OpenFileDialog();
            if (ofd.ShowDialog() == DialogResult.OK)
            {
                Bitmap bmp = new Bitmap(ofd.FileName); // // Загружаем изображение.
                /*Color[,] color = new Color[bmp.Width, bmp.Height];
                for (int y = 0; y < bmp.Height; y++)
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        color[x, y] = bmp.GetPixel(x, y);
                    }
                byte[,] R = new byte[bmp.Size.Width, bmp.Size.Height];
                byte[,] G = new byte[bmp.Size.Width, bmp.Size.Height];
                byte[,] B = new byte[bmp.Size.Width, bmp.Size.Height];
                for (int y = 0; y < bmp.Height; y++)
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        R[x, y] = color[x, y].R;
                        G[x, y] = color[x, y].G;
                        B[x, y] = color[x, y].B;

                    }
                for (int i = 0; i < bmp.Height; i++)
                    for (int j = 0; j < bmp.Width; j++)
                    {
                        bmp.SetPixel(j, i, Color.FromArgb((R[j, i] + G[j, i] + B[j, i]) / 3, (R[j, i] + G[j, i] + B[j, i]) / 3, (R[j, i] + G[j, i] + B[j, i]) / 3));
                    }*/
                pictureBox2.Image = bmp; // Устанавливаем в PictureBox/            

            }

        }

        private void pictureBox2_Click(object sender, EventArgs e)
        {

        }
    }
}
