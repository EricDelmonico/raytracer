using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Microsoft.Xna.Framework;
using System.Runtime.InteropServices;
using System.Threading;

namespace raytracer
{
    public delegate Color PixelFunction(int xPos, int yPos);

    public struct HitRecord
    {
        public static HitRecord Empty;

        public Vec3 p;
        public Vec3 normal;
        public Material mat;
        public double t;
        public double u;
        public double v;
        public bool frontFace;

        public void SetFaceNormal(ref DRay r, ref Vec3 outwardNormal)
        {
            frontFace = Vec3.Dot(r.Direction, outwardNormal) < 0;
            normal = frontFace ? outwardNormal : -outwardNormal;
        }
    }

    public abstract class Hittable
    {
        public abstract bool Hit(ref DRay r, double tMin, double tMax, ref HitRecord rec);
        public abstract bool BoundingBox(double t0, double t1, ref AABB outputBox);
    }

    public class Sphere : Hittable
    {
        public Vec3 Center { get { return center; } }
        public double Radius { get { return radius; } }

        private Vec3 center;
        private double radius;
        public Material mat;

        public Sphere(Vec3 center, double radius, Material mat)
        {
            this.center = center;
            this.radius = radius;
            this.mat = mat;
        }

        public Sphere()
        {
            center = new Vec3(0, 0, -1);
            radius = 0.5;
        }

        public override bool Hit(ref DRay r, double tMin, double tMax, ref HitRecord rec)
        {
            Vec3 oc = r.Position - center;
            double a = Vec3.Dot(r.Direction, r.Direction);
            double half_b = Vec3.Dot(oc, r.Direction);
            double c = Vec3.Dot(oc, oc) - radius * radius;
            double discriminant = half_b * half_b - a * c;
            if (discriminant > 0)
            {
                double root = Math.Sqrt(discriminant);
                double tempT = (-half_b - root) / a;
                if (tempT < tMax && tempT > tMin)
                {
                    rec.t = tempT;
                    rec.p = r.At(rec.t);
                    rec.normal = (rec.p - center) / radius;
                    Vec3 outwardNormal = (rec.p - center) / radius;
                    rec.SetFaceNormal(ref r, ref outwardNormal);
                    rec.mat = mat;
                    GetSphereUV((rec.p - center) / radius, ref rec.u, ref rec.v);
                    return true;
                }
                tempT = (-half_b + root) / a;
                if (tempT < tMax && tempT > tMin)
                {
                    rec.t = tempT;
                    rec.p = r.At(rec.t);
                    rec.normal = (rec.p - center) / radius;
                    Vec3 outwardNormal = (rec.p - center) / radius;
                    rec.SetFaceNormal(ref r, ref outwardNormal);
                    rec.mat = mat;
                    GetSphereUV((rec.p - center) / radius, ref rec.u, ref rec.v);
                    return true;
                }
            }
            return false;
        }

        public override bool BoundingBox(double t0, double t1, ref AABB outputBox)
        {
            outputBox = new AABB(center - new Vec3(radius, radius, radius),
                                 center + new Vec3(radius, radius, radius));
            return true;
        }

        public void GetSphereUV(Vec3 p, ref double u, ref double v)
        {
            double phi = Math.Atan2(p.Z, p.X);
            double theta = Math.Asin(p.Y);
            double pi = Math.PI;
            u = 1 - (phi + pi) / (2 * pi);
            v = (theta + pi / 2) / pi;
        }
    }

    public class MovingSphere : Hittable
    {
        private Vec3 center0, center1;
        private double time0, time1;
        private double radius;
        private Material mat;

        public MovingSphere(Vec3 center0, Vec3 center1, double t0, double t1, double r, Material material)
        {
            this.center0 = center0;
            this.center1 = center1;
            this.time0 = t0;
            this.time1 = t1;
            this.radius = r;
            this.mat = material;
        }

        public override bool BoundingBox(double t0, double t1, ref AABB outputBox)
        {
            AABB box0 = new AABB(Center(t0) - new Vec3(radius, radius, radius),
                                 Center(t0) + new Vec3(radius, radius, radius));
            AABB box1 = new AABB(Center(t1) - new Vec3(radius, radius, radius),
                                 Center(t1) + new Vec3(radius, radius, radius));

            outputBox = Raytracing.SurroundingBox(box0, box1);
            return true;
        }

        public Vec3 Center(double time)
        {
            return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
        }

        public override bool Hit(ref DRay r, double tMin, double tMax, ref HitRecord rec)
        {
            Vec3 oc = r.Position - Center(r.Time);
            double a = Vec3.Dot(r.Direction, r.Direction);
            double half_b = Vec3.Dot(oc, r.Direction);
            double c = Vec3.Dot(oc, oc) - radius * radius;
            double discriminant = half_b * half_b - a * c;
            if (discriminant > 0)
            {
                double root = Math.Sqrt(discriminant);
                double tempT = (-half_b - root) / a;
                if (tempT < tMax && tempT > tMin)
                {
                    rec.t = tempT;
                    rec.p = r.At(rec.t);
                    rec.normal = (rec.p - Center(r.Time)) / radius;
                    Vec3 outwardNormal = (rec.p - Center(r.Time)) / radius;
                    rec.SetFaceNormal(ref r, ref outwardNormal);
                    rec.mat = mat;
                    return true;
                }
                tempT = (-half_b + root) / a;
                if (tempT < tMax && tempT > tMin)
                {
                    rec.t = tempT;
                    rec.p = r.At(rec.t);
                    rec.normal = (rec.p - Center(r.Time)) / radius;
                    Vec3 outwardNormal = (rec.p - Center(r.Time)) / radius;
                    rec.SetFaceNormal(ref r, ref outwardNormal);
                    rec.mat = mat;
                    return true;
                }
            }
            return false;
        }
    }

    public class HittableList : Hittable
    {
        public List<Hittable> objects;

        public HittableList()
        {
            objects = new List<Hittable>();
        }

        public void Clear()
        {
            objects.Clear();
        }

        public void AddRef(ref Hittable obj)
        {
            objects.Add(obj);
        }
        public void Add(Hittable obj)
        {
            objects.Add(obj);
        }

        public override bool Hit(ref DRay r, double tMin, double tMax, ref HitRecord rec)
        {
            HitRecord tempRecord = HitRecord.Empty;
            bool hitAnything = false;
            double closestSoFar = tMax;

            foreach (Hittable obj in objects)
            {
                if (obj.Hit(ref r, tMin, closestSoFar, ref tempRecord))
                {
                    hitAnything = true;
                    closestSoFar = tempRecord.t;
                    rec = tempRecord;
                }
            }

            return hitAnything;
        }

        public override bool BoundingBox(double t0, double t1, ref AABB outputBox)
        {
            if (objects.Count == 0) return false;

            AABB temp = null;
            bool firstBox = true;

            foreach (Hittable obj in objects)
            {
                if (!obj.BoundingBox(t0, t1, ref temp)) return false;
                outputBox = firstBox ? temp : Raytracing.SurroundingBox(outputBox, temp);
                firstBox = false;
            }

            return true;
        }
    }

    public abstract class Material
    {
        public abstract bool Scatter(ref DRay rayIn, ref HitRecord rec, ref Vec3 attenuation, ref DRay scattered);
    }

    public abstract class Texture
    {
        public abstract Vec3 Value(double u, double v, Vec3 p);
    }

    public struct Vec3
    {
        public double this[int i]
        {
            get
            {
                return (new double[] { this.X, this.Y, this.Z})[i];
            }
        }

        public double X;
        public double Y;
        public double Z;

        public static Vec3 Zero
        {
            get
            {
                return new Vec3(0);
            }
        }
        public Vec3 normalized
        {
            get
            {
                return this / Length();
            }
        }

        public Vec3(double number)
        {
            X = number;
            Y = number;
            Z = number;
        }

        public Vec3(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        // functions that operate on this vec specifically
        public void Normalize()
        {
            this /= Length();
        }

        public double LengthSquared()
        {
            return (X * X) + (Y * Y) + (Z * Z);
        }

        public double Length()
        {
            return Math.Sqrt(LengthSquared());
        }

        // helper equations
        public static double Dot(Vec3 a, Vec3 b)
        {
            return (a.X * b.X) + (a.Y * b.Y) + (a.Z * b.Z);
        }

        public static Vec3 Cross(Vec3 a, Vec3 b)
        {
            return new Vec3(a.Y * b.Z - a.Z * b.Y,
                            a.Z * b.X - a.X * b.Z,
                            a.X * b.Y - a.Y * b.X);
        }

        // operators
        public static Vec3 operator -(Vec3 a)
        {
            return new Vec3(-a.X, -a.Y, -a.Z);
        }

        public static Vec3 operator +(Vec3 a, Vec3 b)
        {
            return new Vec3(a.X + b.X,
                            a.Y + b.Y,
                            a.Z + b.Z);
        }

        public static Vec3 operator -(Vec3 a, Vec3 b)
        {
            return new Vec3(a.X - b.X,
                            a.Y - b.Y,
                            a.Z - b.Z);
        }

        public static Vec3 operator *(Vec3 a, double t)
        {
            return new Vec3(a.X * t, a.Y * t, a.Z * t);
        }

        public static Vec3 operator *(double t, Vec3 a)
        {
            return new Vec3(a.X * t, a.Y * t, a.Z * t);
        }

        public static Vec3 operator /(Vec3 a, double t)
        {
            return a * (1 / t);
        }

        public static Vec3 operator *(Vec3 a, Vec3 b)
        {
            return new Vec3((a.X * b.X), (a.Y * b.Y), (a.Z * b.Z));
        }
    }

    public struct DRay
    {
        public Vec3 Position { get { return position; } }
        public Vec3 Direction { get { return direction; } }
        public double Time { get { return time; } }

        private Vec3 position;
        private Vec3 direction;
        private double time;

        public DRay(Vec3 position, Vec3 direction, double time = 0.0)
        {
            this.position = position;
            this.direction = direction;
            this.time = time;
        }

        // p = p + vt
        public Vec3 At(double time)
        {
            return position + (direction * time);
        }
    }

    public class AABB
    {
        public Vec3 Min { get { return min; } }
        public Vec3 Max { get { return max; } }

        private Vec3 min;
        private Vec3 max;

        public AABB(Vec3 min, Vec3 max)
        {
            this.min = min;
            this.max = max;
        }

        public bool Hit(ref DRay r, double tmin, double tmax)
        {
            for (int a = 0; a < 3; a++)
            {
                double invD = 1.0f / r.Direction[a];
                double t0 = (min[a] - r.Position[a]) * invD;
                double t1 = (max[a] - r.Position[a]) * invD;
                if (invD < 0.0f)
                    Raytracing.QuickSwap(ref t0, ref t1);
                tmin = t0 > tmin ? t0 : tmin;
                tmax = t1 < tmax ? t1 : tmax;
                if (tmax <= tmin)
                    return false;
            }
            return true;
        }
    }

    public class BVHNode : Hittable
    {
        private delegate int BVHCompare(Hittable a, Hittable b);
        private class BVHComparer : IComparer<Hittable>
        {
            public BVHCompare CompareFunc { set { Comparer = value; } }

            private BVHCompare Comparer;

            public int Compare(Hittable x, Hittable y)
            {
                return Comparer(x, y);
            }
        }

        public Hittable Right { get { return right; } }
        public Hittable Left { get { return left; } }
        public AABB Box { get { return box; } }

        private Hittable left;
        private Hittable right;
        private AABB box;

        public BVHNode(HittableList list, double time0, double time1) : this(list.objects, 0, list.objects.Count, time0, time1) { }
        public BVHNode(List<Hittable> objects, int start, int end, double time0, double time1)
        {
            int axis = (int)Raytracing.RandomInRange(0, 2);
            BVHComparer comparator = new BVHComparer();
            switch (axis)
            {
                case 0:
                    comparator.CompareFunc = BoxCompareX;
                    break;
                case 1:
                    comparator.CompareFunc = BoxCompareY;
                    break;
                case 2:
                    comparator.CompareFunc = BoxCompareZ;
                    break;
                default:
                    break;
            }

            int objectSpan = end - start;

            if (objectSpan == 1)
            {
                left = right = objects[start];
            }
            else if (objectSpan == 2)
            {
                if (comparator.Compare(objects[start], objects[start + 1]) > 0)
                {
                    left = objects[start];
                    right = objects[start + 1];
                }
                else
                {
                    left = objects[start + 1];
                    right = objects[start];
                }
            }
            else
            {
                objects.Sort(start, end - start, comparator);

                int mid = start + objectSpan / 2;
                left = new BVHNode(objects, start, mid, time0, time1);
                right = new BVHNode(objects, mid, end, time0, time1);
            }

            AABB boxLeft = null, boxRight = null;

            if (!left.BoundingBox(time0, time1, ref boxLeft)
               || !right.BoundingBox(time0, time1, ref boxRight))
            {
                throw new ArgumentNullException("no bounding box in bvh node constructor");
            }

            box = Raytracing.SurroundingBox(boxLeft, boxRight);
        }

        public override bool BoundingBox(double t0, double t1, ref AABB outputBox)
        {
            outputBox = box;
            return true;
        }

        public override bool Hit(ref DRay r, double tMin, double tMax, ref HitRecord rec)
        {
            if (!box.Hit(ref r, tMin, tMax))
            {
                return false;
            }

            bool hitL = left.Hit(ref r, tMin, tMax, ref rec);
            bool hitR = right.Hit(ref r, tMin, tMax, ref rec);

            return hitL || hitR;
        }

        private int BoxCompare(Hittable a, Hittable b, int axis)
        {
            AABB boxA = null, boxB = null;

            if (!a.BoundingBox(0, 0, ref boxA) || !b.BoundingBox(0, 0, ref boxB))
            {
                throw new ArgumentNullException("no bounding box in bvh node constructor");
            }

            return boxA.Min[axis] < boxB.Min[axis] ? 0 : 1;
        }
        private int BoxCompareX(Hittable a, Hittable b)
        {
            return BoxCompare(a, b, 0);
        }
        private int BoxCompareY(Hittable a, Hittable b)
        {
            return BoxCompare(a, b, 1);
        }
        private int BoxCompareZ(Hittable a, Hittable b)
        {
            return BoxCompare(a, b, 1);
        }
    }

    static class Raytracing
    {
        #region materials
        public class Lambertian : Material
        {
            private Texture albedo;

            public Lambertian(Texture albedo)
            {
                this.albedo = albedo;
            }

            public override bool Scatter(ref DRay rayIn, ref HitRecord rec, ref Vec3 attenuation, ref DRay scattered)
            {
                // more accurate technically: use RandomInHemisphere(ref rec.normal)
                Vec3 scatterDirection = rec.normal + RandomInUnitSphere();
                scattered = new DRay(rec.p, scatterDirection, rayIn.Time);
                attenuation = albedo.Value(rec.u, rec.v, rec.p);
                return true;
            }
        }
        public class Metal : Material
        {
            private Vec3 albedo;
            private double fuzziness;

            public Metal(Color albedo, double fuzziness)
            {
                this.albedo = new Vec3(albedo.R / 255.0, albedo.G / 255.0, albedo.B / 255.0);
                this.fuzziness = fuzziness;
            }

            public override bool Scatter(ref DRay rayIn, ref HitRecord rec, ref Vec3 attenuation, ref DRay scattered)
            {
                Vec3 direction = rayIn.Direction;
                direction.Normalize();
                Vec3 reflected = Reflect(ref direction, ref rec.normal);
                scattered = new DRay(rec.p, reflected + fuzziness * RandomInUnitSphere());
                attenuation = albedo;
                return Vec3.Dot(scattered.Direction, rec.normal) > 0;
            }
        }
        public class Dielectric : Material
        {
            public double IndexOfRefraction { get { return refInd; } }

            private double refInd;

            public Dielectric(double rIndex)
            {
                refInd = rIndex;
            }

            public override bool Scatter(ref DRay rayIn, ref HitRecord rec, ref Vec3 attenuation, ref DRay scattered)
            {
                attenuation = new Vec3(1.0, 1.0, 1.0);
                double etaOverEtaPrime;
                if (rec.frontFace)
                {
                    etaOverEtaPrime = 1 / refInd;
                }
                else
                {
                    etaOverEtaPrime = refInd;
                }

                Vec3 unitDirection = rayIn.Direction.normalized;
                double cosTheta = Math.Min(Vec3.Dot(-unitDirection, rec.normal), 1.0);
                double sinTheta = Math.Sqrt(1.0 - cosTheta * cosTheta);
                if (etaOverEtaPrime * sinTheta > 1.0)
                {
                    Vec3 reflected = Reflect(ref unitDirection, ref rec.normal);
                    scattered = new DRay(rec.p, reflected);
                    return true;
                }
                double reflectProb = Schlick(cosTheta, etaOverEtaPrime);
                if (NextRandom() < reflectProb)
                {
                    Vec3 reflected = Reflect(ref unitDirection, ref rec.normal);
                    scattered = new DRay(rec.p, reflected);
                    return true;
                }

                Vec3 refracted = Refract(ref unitDirection, ref rec.normal, etaOverEtaPrime);
                scattered = new DRay(rec.p, refracted);
                return true;
            }
        }
        #endregion

        #region textures

        public class SolidColor : Texture
        {
            private Vec3 color;

            public SolidColor(Vec3 color)
            {
                this.color = color;
            }
            public SolidColor(double red, double green, double blue) : this(new Vec3(red, green, blue)) { }

            public override Vec3 Value(double u, double v, Vec3 p)
            {
                return color;
            }
        }

        public class CheckerTexture : Texture
        {
            private Texture even;
            private Texture odd;

            public CheckerTexture(Texture t0, Texture t1)
            {
                even = t0;
                odd = t1;
            }

            public override Vec3 Value(double u, double v, Vec3 p)
            {
                double sines = Math.Sin(10 * p.X) * Math.Sin(10 * p.Y) * Math.Sin(10 * p.Z);
                if (sines < 0)
                {
                    return odd.Value(u, v, p);
                }
                else
                {
                    return even.Value(u, v, p);
                }
            }
        }

        #endregion

        #region camera stuff
        private const double deg2Rad = Math.PI / 180.0;
        private const double ASPECT_RATIO = 16.0 / 9.0;

        public class Camera
        {
            Vec3 origin;
            Vec3 horizontal;
            Vec3 vertical;
            Vec3 lowerLeftCorner;

            double lensRadius;
            Vec3 u;
            Vec3 v;
            Vec3 w;

            double time0;
            double time1;

            public Camera(Vec3 lookfrom,
                          Vec3 lookat,
                          Vec3 up,
                          double vertFOV,
                          double aspectRatio,
                          double aperture,
                          double focusDistance,
                          double t0 = 0,
                          double t1 = 0)
            {
                double theta = vertFOV * deg2Rad;
                double h = Math.Tan(theta / 2);
                double viewportHeight = 2.0 * h;
                double viewportWidth = aspectRatio * viewportHeight;

                w = (lookfrom - lookat).normalized;
                u = Vec3.Cross(up, w).normalized;
                v = Vec3.Cross(w, u);

                origin = lookfrom;
                horizontal = focusDistance * viewportWidth * u;
                vertical = focusDistance * viewportHeight * v;
                lowerLeftCorner = origin - horizontal / 2 - vertical / 2 - focusDistance * w;

                lensRadius = aperture / 2;

                this.time0 = t0;
                this.time1 = t1;
            }

            public DRay GetRay(double horizontalAmount, double verticalAmount)
            {
                Vec3 rd = lensRadius * RandomInUnitDisk();
                Vec3 offset = u * rd.X + v * rd.Y;

                return new DRay(origin + offset,
                                lowerLeftCorner + horizontalAmount * horizontal + verticalAmount * vertical - origin - offset,
                                RandomInRange(time0, time1));
            }
        }

        #endregion

        #region raytracing members
        // dimensions of the final image
        public static int IMG_WIDTH = 256;
        public static int IMG_HEIGHT = (int)(IMG_WIDTH / ASPECT_RATIO);
        public static Random rng;

        public static int scanlinesRemaining = IMG_HEIGHT;
        public static int maxThreads = 4;

        // the function that will operate on every pixel of the final image
        private static readonly PixelFunction pixelFunction = Raytrace;

        // the world full of hittable objects
        private static HittableList world;
        private static HittableList randomScene;
        private static Camera camera;
        private static double samplesPerPixel;

        // the max depth of light bouces
        private const int maxDepth = 50;

        public static Color[] data;

        private static double[] randomValues;
        private static int currentRandom;
        private const int maxRandom = 1000000;
        private const int error = 1000;

        private static HittableList SCENE;
        #endregion

        public static void Init()
        {
            rng = new Random();
            // populating an array of random values
            randomValues = new double[maxRandom + error];
            for (int i = 0; i < maxRandom + error; i++)
            {
                randomValues[i] = rng.NextDouble();
            }

            #region creating world hittablelist
            world = new HittableList();
            world.Add(new Sphere(new Vec3(0, 0, -1), 0.5, new Lambertian(new SolidColor(0.1, 0.2, 0.5))));
            world.Add(new Sphere(new Vec3(0, -100.5, -1), 100, new Lambertian(new SolidColor(0.8, 0.8, 0.0))));
            world.Add(new Sphere(new Vec3(-1, 0, -1), 0.5, new Dielectric(1.5)));
            world.Add(new Sphere(new Vec3(-1, 0, -1), -0.45, new Dielectric(1.5)));
            world.Add(new Sphere(new Vec3(1, 0, -1), 0.5, new Metal(new Color(.8f, .6f, .2f), 0)));

            // a red and a blue sphere, for testing
            //double R = Math.Cos(Math.PI / 4);
            //world.Add(new Sphere(new Vec3(-R, 0, -1), R, new Lambertian(new Color(0, 0, 1f))));
            //world.Add(new Sphere(new Vec3(R, 0, -1), R, new Lambertian(new Color(1f, 0, 0))));
            #endregion

            #region creating randomScene hittablelist
            randomScene = new HittableList();
            Texture checker = new CheckerTexture(new SolidColor(0.2, 0.3, 0.1),
                                                 new SolidColor(0.9, 0.9, 0.9));
            Material groundMat = new Lambertian(checker);
            randomScene.Add(new Sphere(new Vec3(0, -1000, 0), 1000, groundMat));

            double ssr = 0.2; // small spheres' radius

            for (int i = -11; i < 11; i++)
            {
                for (int j = -11; j < 11; j++)
                {
                    double chooseMat = NextRandom();
                    Vec3 center = new Vec3(i + 0.9 * NextRandom(), 0.2, j + 0.9 * NextRandom());

                    if ((center - new Vec3(4, 0.2, 0)).Length() > 0.9)
                    {
                        Material sphereMat;

                        if (chooseMat < 0.8)
                        {
                            //diffuse
                            Vec3 albedo = RandomVec() * RandomVec();
                            sphereMat = new Lambertian(new SolidColor(albedo));
                            Vec3 center2 = center + new Vec3(0, RandomInRange(0.0, 0.5), 0);
                            randomScene.Add(new MovingSphere(center, center2, 0, 1, ssr, sphereMat));
                        }
                        else if (chooseMat < 0.95)
                        {
                            //metal
                            Color albedo = Vec3ToColor(RandomVec(0.5, 1));
                            double fuzz = RandomInRange(0, 0.5);
                            sphereMat = new Metal(albedo, fuzz);
                            randomScene.Add(new Sphere(center, ssr, sphereMat));
                        }
                        else
                        {
                            //glass
                            sphereMat = new Dielectric(1.5);
                            randomScene.Add(new Sphere(center, ssr, sphereMat));
                        }
                    }
                }
            }

            Material material1 = new Dielectric(1.5);
            randomScene.Add(new Sphere(new Vec3(0, 1, 0), 1.0, material1));

            Material material2 = new Lambertian(new SolidColor(0.4, 0.2, 0.1));
            randomScene.Add(new Sphere(new Vec3(-4, 1, 0), 1.0, material2));

            Material material3 = new Metal(new Color(0.7f, 0.6f, 0.5f), 0.0);
            randomScene.Add(new Sphere(new Vec3(4, 1, 0), 1.0, material3));

            HittableList randomSceneBVH = new HittableList();
            randomSceneBVH.Add(new BVHNode(randomScene, 0, 1));
            #endregion

            #region creating camera
            Vec3 lookfrom = new Vec3(13, 2, 3);
            Vec3 lookat = new Vec3(0, 0, 0);
            Vec3 up = new Vec3(0, 1, 0);
            double fov = 20;
            double aspectRatio = ASPECT_RATIO;
            double aperture = 0.1;
            double focusDistance = 10.0;

            camera = new Camera(lookfrom,
                                lookat,
                                up,
                                fov,
                                aspectRatio,
                                aperture,
                                focusDistance,
                                0,
                                1);
            #endregion

            samplesPerPixel = 100;
            data = new Color[IMG_WIDTH * IMG_HEIGHT];
            SCENE = randomSceneBVH;
        }

        // ------------------------
        // MAIN RENDERING FUNCTIONS
        // ------------------------

        /// <summary>
        /// Gets the data for the final render using the pixelFunction, and returns it
        /// </summary>
        /// <returns>Color data for use in Texture2D.SetData function</returns>
        public static Color[] GetTextureData()
        {
            // single threaded version
            for (int i = IMG_HEIGHT; i > 0; i--)
            {
                for (int j = 0; j < IMG_WIDTH; j++)
                {
                    data[(IMG_HEIGHT - i) * IMG_WIDTH + j] = pixelFunction(j, i);
                }
                scanlinesRemaining--;
            }

            return data;
        }

        public static Color[] GetTextureDataInSpiral()
        {
            int err = 1;
            int maxX = IMG_WIDTH - err, maxY = IMG_HEIGHT - err, minX = 0 + err, minY = 0 + err;
            int cX = 1 + err, cY = 1 + err; // current x and current y
            int cD = 1; // current direction: NESW = 1234
            while (maxX > minX || maxY > minY)
            {
                if (cX < maxX && cD == 1)
                {
                    cX += err;
                    data[(IMG_HEIGHT - cY) * IMG_WIDTH + cX] = pixelFunction(cX, cY);
                    if (cX >= maxX)
                    {
                        maxX -= err;
                        cD = 2;
                    }
                    continue;
                }
                if (cY > minY && cD == 2)
                {
                    cY -= err;
                    data[(IMG_HEIGHT - cY) * IMG_WIDTH + cX] = pixelFunction(cX, cY);
                    if (cY <= minY)
                    {
                        minY += err;
                        cD = 3;
                    }
                    continue;
                }
                if (cX > minX && cD == 3)
                {
                    cX -= err;
                    data[(IMG_HEIGHT - cY) * IMG_WIDTH + cX] = pixelFunction(cX, cY);
                    if (cX <= minX)
                    {
                        minX += err;
                        cD = 4;
                    }
                    continue;
                }
                if (cY < maxY && cD == 4)
                {
                    cY += err;
                    data[(IMG_HEIGHT - cY) * IMG_WIDTH + cX] = pixelFunction(cX, cY);
                    if (cY >= maxY)
                    {
                        maxY -= err;
                        cD = 1;
                    }
                    continue;
                }
            }

            return data;
        }

        public static Color[] GetTextureDataThreadedByBoxes()
        {
            // making the event that will be waited for
            // in order to make sure all threads finish
            // before pushing the texture data to Game1
            int threadCount = 1024;
            CountdownEvent countdownEvent = new CountdownEvent(threadCount);
            // rendering threadCount boxes and returning the data
            for (int box = 0; box < threadCount; box++)
            {
                // splitting the render into threadCount boxes,
                // each of which will be sent to the thread pool
                int width = IMG_WIDTH / 32;
                int height = IMG_HEIGHT / 32;
                int startX = (box % 32) * width;
                int startY = (box / 32) * height;

                ThreadPool.QueueUserWorkItem(
                    (Object stateInfo) =>
                    {
                        for (int i = height + startY; i > startY; i--)
                        {
                            for (int j = startX; j - startX < width; j++)
                            {
                                //// checkerboard ish? fun to play with but dont keep lol
                                //if (j*i % 5 != 0)
                                //{
                                //    data[(IMG_HEIGHT - i) * IMG_WIDTH + j] = data[(IMG_HEIGHT - i) * IMG_WIDTH + j - 1];
                                //    continue;
                                //}
                                data[(IMG_HEIGHT - i) * IMG_WIDTH + j] = pixelFunction(j, i);
                            }
                        }
                        countdownEvent.Signal();
                    });
            }
            countdownEvent.Wait();

            return data;
        }

        public static Color[] GetTextureDataThreadedByScanline()
        {
            for (int i = IMG_HEIGHT; i > 0; i--)
            {
                if (maxThreads <= 0)
                {
                    i++;
                    continue;
                }
                ThreadPool.QueueUserWorkItem(new WaitCallback(
                    (object stateInfo) =>
                    {
                        for (int j = 0; j < IMG_WIDTH; j++)
                        {
                            data[(IMG_HEIGHT - (int)stateInfo) * IMG_WIDTH + j] = pixelFunction(j, (int)stateInfo);
                        }
                        scanlinesRemaining--;
                        maxThreads++;
                    }
                    ), i);
                maxThreads--;
            }
            return data;
        }

        public static Color Raytrace(int x, int y)
        {
            // set this pixel's color
            Vec3 pixelColor = Vec3.Zero;
            for (int s = 0; s < samplesPerPixel; s++)
            {
                double rX = NextRandom();
                double rY = NextRandom();
                double horizontalAmount = (x + rX) / (IMG_WIDTH - 1);
                double verticalAmount = (y + rY) / (IMG_HEIGHT - 1);
                DRay r = camera.GetRay(horizontalAmount, verticalAmount);
                pixelColor += RayColor(ref r, ref SCENE, maxDepth);
            }

            double scale = 1.0 / samplesPerPixel;
            pixelColor.X = Math.Sqrt(scale * pixelColor.X);
            pixelColor.Y = Math.Sqrt(scale * pixelColor.Y);
            pixelColor.Z = Math.Sqrt(scale * pixelColor.Z);
            pixelColor.X = Clamp(pixelColor.X, 0, 0.999);
            pixelColor.Y = Clamp(pixelColor.Y, 0, 0.999);
            pixelColor.Z = Clamp(pixelColor.Z, 0, 0.999);

            return Vec3ToColor(pixelColor);
        }

        // ----------------
        // HELPER FUNCTIONS
        // ----------------

        public static AABB SurroundingBox(AABB box0, AABB box1)
        {
            // constructing a box that surrounds both the inputted boxes
            Vec3 min = new Vec3(Math.Min(box0.Min.X, box1.Min.X),
                                Math.Min(box0.Min.Y, box1.Min.Y),
                                Math.Min(box0.Min.Z, box1.Min.Z));

            Vec3 max = new Vec3(Math.Max(box0.Max.X, box1.Max.X),
                                Math.Max(box0.Max.Y, box1.Max.Y),
                                Math.Max(box0.Max.Z, box1.Max.Z));

            return new AABB(min, max);
        }

        public static void QuickSwap(ref double x, ref double y)
        {
            double t = x;
            x = y;
            y = t;
        }

        private static double Clamp(double value, double min, double max)
        {
            if (value < min)
                return min;
            if (value > max)
                return max;

            return value;
        }
        private static Color Vec3ToColor(Vec3 color)
        {
            return new Color((float)color.X, (float)color.Y, (float)color.Z);
        }

        private static double Schlick(double cos, double refInd)
        {
            double r0 = (1 - refInd) / (1 + refInd);
            r0 = r0 * r0;
            return r0 + (1 - r0) * Math.Pow((1 - cos), 5);
        }

        private static Vec3 Reflect(ref Vec3 v, ref Vec3 n)
        {
            // v pointing inwards, with "-" pointing the reflected ray outwards,
            // v dot n being the component of v in the direction of n, and multiplying
            // by n to get the vector that, when multiplied by 2, will yield reflection
            // \  n  /|
            //  \ | / | <-- b on the outside, reflected v on the inside
            //v->\|/  |
            //--------------
            // also\  | 
            //  v-> \ | <-- b
            //       \| 

            return v - 2 * Vec3.Dot(v, n) * n;
        }

        private static Vec3 Refract(ref Vec3 uv, ref Vec3 n, double etaOverEtaPrime)
        {
            double cosTheta = Vec3.Dot(-uv, n);
            // component of r parallel with the normal
            Vec3 rOutParallel = etaOverEtaPrime * (uv + cosTheta * n);
            // component of r perpindicular with the normal
            Vec3 rOutPerp = -Math.Sqrt(1.0 - rOutParallel.LengthSquared()) * n;
            return rOutParallel + rOutPerp;
        }

        private static Vec3 RandomInUnitDisk()
        {
            while (true)
            {
                Vec3 p = new Vec3(RandomInRange(-1, 1), RandomInRange(-1, 1), 0);
                if (p.LengthSquared() >= 1)
                {
                    continue;
                }

                return p;
            }
        }

        private static Vec3 RandomVec(double min = 0, double max = 1.0)
        {
            double range = max - min;
            double rand1 = NextRandom();
            double rand2 = NextRandom();
            double rand3 = NextRandom();
            double x = min + (rand1 * range);
            double y = min + (rand2 * range);
            double z = min + (rand3 * range);
            return new Vec3(x, y, z);
        }

        private static Vec3 RandomInHemisphere(ref Vec3 normal)
        {
            Vec3 inUnitSphere = RandomInUnitSphere();
            if (Vec3.Dot(inUnitSphere, normal) > 0.0)
            {
                return inUnitSphere;
            }
            else
            {
                return -inUnitSphere;
            }
        }

        private static Vec3 RandomUnitVector()
        {
            double a = RandomInRange(0, Math.PI * 2);
            double z = RandomInRange(-1, 1);
            double r = Math.Sqrt(1 - z * z);
            return new Vec3((r * Math.Cos(a)), (r * Math.Sin(a)), (z));
        }

        private static Vec3 RandomInUnitSphere()
        {
            while (true)
            {
                Vec3 p = RandomVec(-1.0, 1.0);
                if (p.LengthSquared() >= 1)
                {
                    continue;
                }

                return p;
            }
        }

        public static double RandomInRange(double min, double max)
        {
            return (min + (NextRandom() * (max - min)));
        }

        /// <summary>
        /// creates a simple blue->white background gradient in the
        /// case of rays that don't intersect anything in the scene
        /// </summary>
        /// <param name="r">the ray to test a color for</param>
        /// <returns>the color a background pixel should be according to the ray's direction</returns>
        private static Vec3 RayColor(ref DRay r, ref HittableList world, int depth)
        {
            if (depth <= 0)
            {
                return Vec3.Zero;
            }

            HitRecord rec = HitRecord.Empty;

            // get the ray's normalized direction
            Vec3 dir = r.Direction.normalized;
            if (world.Hit(ref r, 0.00001, double.MaxValue, ref rec))
            {
                // material method
                DRay scattered = new DRay();
                Vec3 attenuation = Vec3.Zero;
                if (rec.mat.Scatter(ref r, ref rec, ref attenuation, ref scattered))
                {
                    return attenuation * RayColor(ref scattered, ref world, depth - 1);
                }

                // no material method: 
                // hemisphere diffuse: Vec3 target = rec.p + RandomInHemisphere(ref rec.normal);
                // normal diffuse
                //Vec3 target = rec.p + rec.normal + RandomUnitVector();
                //DRay rayContinue = new DRay(rec.p, target - rec.p);
                //return (0.5f * RayColor(ref rayContinue, ref world, depth - 1));

                // normal-based color gradient
                //return (0.5f * (rec.normal + new Vec3(1.0f, 1.0f, 1.0f)));
            }

            // convert the y-dir to a value from 0-1
            double y = (dir.Y + 1.0) / 2.0;

            // lerping white value to off-blue to make the gradient
            Vec3 white = new Vec3(1.0, 1.0, 1.0);
            Vec3 offBlue = new Vec3(0.5, 0.7, 1.0);

            return white * (1.0 - y) + offBlue * y;
        }

        // -----------------
        // TESTING FUNCTIONS
        // -----------------

        /// <summary>
        /// Generates a simple gradient image to display on the 
        /// screen to check if rendering is working properly
        /// </summary>
        /// <param name="pixelX">the x-coordinate of the pixel to operate on</param>
        /// <param name="pixelY">the y-coordinate of the pixel to operate on</param>
        /// <returns>A Color in the gradient</returns>
        private static Color GenerateSampleImage(int pixelX, int pixelY)
        {
            Color pixelColor = new Color((pixelX / IMG_WIDTH),
                                         (pixelY / IMG_HEIGHT),
                                         0.25f);

            return pixelColor;
        }

        private static double NextRandom()
        {
            if (currentRandom >= maxRandom)
            {
                currentRandom = 0;
            }
            return randomValues[currentRandom++];
        }
    }
}
