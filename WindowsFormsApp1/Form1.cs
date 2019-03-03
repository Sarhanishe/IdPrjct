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
        double[] forest1 = { 0.7, 0.94, 0.36, 0.52 };
        double[] city2 = { -25.0088, -23, 3300, 839.4985, 933.3491 };
        double[] river3 = { -2.2365, -1.0706, 29.2263, 65.0149 };
        public CR()
        {
            InitializeComponent();
        }
        //Рекурсия
        public void Identification(double[,] matr, int k)
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
            if (n < 6)
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
        public void Identification(double[,] matr)
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
        public double[] Kf(double[,] pix, int k)
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
            double[,] m = new double[p/s, p/s];
            for (int i = (k / s) * p / s , ii=0; i < ((k / s) + 1) * p / s; i++, ii++)
                for (int j = (k % s) * p / s , jj=0; j < ((k % s) + 1) * p / s; j++ , jj++)
                    m[ii,jj]= pix[i, j] ;
            double[] arr = new double[2];
            var mtr = Matrix<double>.Build.DenseOfArray(m);
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
            arr[0] = ((sing.Length) * sumxy - (sumx * sumy)) / ((sing.Length) * sumx1 - sumx * sumx);
            arr[1] = (sumy - arr[0] * sumx) / ((sing.Length));
            
            if (k == 61 && n == 5)
            {
                using (StreamWriter svdfile = new StreamWriter(@"ftest.txt", append: true))
                {
                    // svdfile.Write(i);
                    svdfile.Write("a0 = ");
                    svdfile.Write(arr[0]);
                    svdfile.Write(" b0 = ");
                    svdfile.Write(arr[1]);
                    svdfile.Write(" n = ");
                    svdfile.Write(p/s);
                    svdfile.WriteLine();
                   /* for (int i = 0; i < p ; i++)
                    {
                        for (int j = 0; j < p ; j++)
                        {
                            svdfile.Write(" ");
                            svdfile.Write(pix[i, j]);
                        }
                        svdfile.WriteLine();
                    }*/
                };
            }
            if (k == 1 && n == 5)
            {
                using (StreamWriter svdfile = new StreamWriter(@"futest.txt", append: true))
                {
                    // svdfile.Write(i);
                    svdfile.Write("a0 = ");
                    svdfile.Write(arr[0]);
                    svdfile.Write(" b0 = ");
                    svdfile.Write(arr[1]);
                    svdfile.Write(" n = ");
                    svdfile.Write(p / s);
                    svdfile.WriteLine();
                    /* for (int i = 0; i < p ; i++)
                     {
                         for (int j = 0; j < p ; j++)
                         {
                             svdfile.Write(" ");
                             svdfile.Write(pix[i, j]);
                         }
                         svdfile.WriteLine();
                     }*/
                };
            }
            if (k == 2 && n == 5)
            {
                using (StreamWriter svdfile = new StreamWriter(@"futest.txt", append: true))
                {
                    // svdfile.Write(i);
                    svdfile.Write("a0 = ");
                    svdfile.Write(arr[0]);
                    svdfile.Write(" b0 = ");
                    svdfile.Write(arr[1]);
                    svdfile.Write(" n = ");
                    svdfile.Write(p / s);
                    svdfile.WriteLine();
                    /* for (int i = 0; i < p ; i++)
                     {
                         for (int j = 0; j < p ; j++)
                         {
                             svdfile.Write(" ");
                             svdfile.Write(pix[i, j]);
                         }
                         svdfile.WriteLine();
                     }*/
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
        public int Eq(double[] kf)
        {
            int n=0;
           /* if (kf[0] == 0 && kf[1] == 0)
                n = 1;*/
            if ((kf[0] >= forest1[0]) && (kf[0] <= forest1[1]) && (kf[1] >= forest1[2]) && (kf[1] <= forest1[3]))
                n = 1;
            else if ((kf[0] >= city2[0]) && (kf[0] <= city2[1]) && (kf[1] >= city2[2]) && (kf[1] <= city2[3]))
                n=2;
            else if ((kf[0] >= river3[0]) && (kf[0] <= river3[1]) && (kf[1] >= river3[2]) && (kf[1] <= river3[3]))
                n=3;
            return n;

        }
        //Закрашивание
        public void Pnt(double[,] pix, int k,int c)
        {
            int p = pix.GetLength(1);
            int n = 1;
            while ((k - System.Math.Pow(4, n)) >= 0)
            {
                k = k - (int)System.Math.Pow(4, n);
                n++;
            }
            int s = (int)System.Math.Pow(2, n);
            for (int i = (k / s) * p / s; i < ((k / s) + 1) * p / s; i++)
                for (int j = (k % s) * p / s; j < ((k % s) + 1) * p / s; j++)
                    pix[i,j] = c;


        }
        //Основная часть
        private void button1_Click(object sender, EventArgs e)
        {
            Bitmap bmp = (Bitmap)pictureBox2.Image;//создаём объект битмапа для записи в него чёрнобелого изображения
            int M = bmp.Height;
            int N = bmp.Width;
            double[,] a = new double[M, N];
            for (int y = 0; y < N; y++)
                for (int x = 0; x < M; x++)
                {
                    a[x, y] = bmp.GetPixel(y, x).R;
                }

           /* double[] koef = new double[2];
            koef = Kf(a);//получение коэффициентов образа

            using (StreamWriter svdfile = new StreamWriter(@"ftest.txt", append: true))
            {
                // svdfile.Write(i);
                svdfile.Write("a0 = ");
                svdfile.Write(koef[0]);
                svdfile.Write(" b0 = ");
                svdfile.Write("lugkjgjhgjhgmjgmjghjjgjh");
                svdfile.WriteLine();
            };*/
            textBox1.Text = (0 % 2).ToString();
            Identification(a);
            textBox1.Text = "end id";
            Bitmap image1 = new Bitmap(bmp.Height, bmp.Width);
            Color newColor;
            // Loop through the images pixels to reset color.
            for (int x = 0; x < image1.Width; x++)
            {
                for (int y = 0; y < image1.Height; y++)
                {
                    switch (a[x,y])
                    {
                        case 1:
                             newColor = Color.FromArgb(0, 250, 0);
                            break;
                        case 2:
                             newColor = Color.FromArgb(250, 0, 0);
                            break;
                        case 3:
                             newColor = Color.FromArgb(0, 0,250);
                            break;
                        default:
                            newColor = Color.FromArgb(0, 0, 0);
                            break;
                    }
                    
                    image1.SetPixel(x, y, newColor);
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

    }
}
