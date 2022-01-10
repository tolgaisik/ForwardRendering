#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"
#include "Matrix4.h"

using namespace tinyxml2;
using namespace std;
void write_ppm(const char *filename, unsigned char *data, int width, int height);

Matrix4 *getViewPortMatrix(Camera *camera)
{
    double viewPort[4][4] = {
        {static_cast<double>(camera->horRes) / 2, 0, 0, (static_cast<double>(camera->horRes) - 1) / 2}, {0, static_cast<double>(camera->verRes) / 2, 0, (static_cast<double>(camera->verRes) - 1) / 2}, {0, 0, 0.5, 0.5}, {0, 0, 0, 1}};
    return new Matrix4(viewPort);
}

Matrix4 *orthographicMatrix(Camera *camera)
{
    double orthographic[4][4] = {
        {2 / (camera->right - camera->left), 0, 0, -(camera->right + camera->left) / (camera->right - camera->left)},
        {0, 2 / (camera->top - camera->bottom), 0, -(camera->top + camera->bottom) / (camera->top - camera->bottom)},
        {0, 0, -2 / (camera->far - camera->near), -(camera->far + camera->near) / (camera->far - camera->near)},
        {0, 0, 0, 1}};
    return new Matrix4(orthographic);
}

Matrix4 *perspectiveMatrix(Camera *camera)
{
    double perspective[4][4] = {{2 * camera->near / (camera->right - camera->left), 0, (camera->right + camera->left) / (camera->right - camera->left), 0},
                                {0, 2 * camera->near / (camera->top - camera->bottom), (camera->top + camera->bottom) / (camera->top - camera->bottom), 0},
                                {0, 0, -(camera->far + camera->near) / (camera->far - camera->near), -2 * camera->far * camera->near / (camera->far - camera->near)},
                                {0, 0, -1, 0}};
    return new Matrix4(perspective);
}

Matrix4 *getCameraTranslationMatrix(Camera *camera)
{
    double innerVector[4][4] = {
        {1, 0, 0, -(camera->pos.x)},
        {0, 1, 0, -(camera->pos.y)},
        {0, 0, 1, -(camera->pos.z)},
        {0, 0, 0, 1}};
    return new Matrix4(innerVector);
}

Matrix4 *getCameraRotationMatrix(Camera *camera)
{
    double rotationArray[4][4] = {
        {camera->u.x, camera->u.y, camera->u.z, 0}, {camera->v.x, camera->v.y, camera->v.z, 0}, {camera->w.x, camera->w.y, camera->w.z, 0}, {0, 0, 0, 1}};
    return new Matrix4(rotationArray);
}

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
    perspective projection, culling, clipping, viewport
*/
// perspective projection, viewport transformation function
Matrix4 *projectionTypeSwitch(Camera *camera)
{
    // perspective projection
    if (camera->projectionType == 0)
    {
        return perspectiveMatrix(camera);
    }
    // orthographic projection
    else if (camera->projectionType == 1)
    {
        return orthographicMatrix(camera);
    }
    // no projection for security reason
    else
    {
        return nullptr;
    }
}

void rotateMesh(Matrix4 &matrix, const Rotation *rotation)
{
    Vec3 u = normalizeVec3(Vec3(rotation->ux, rotation->uy, rotation->uz, 1)), v;
    double smallest_u = std::min(abs(u.x), std::min(abs(u.y), abs(u.z)));
    if (abs(u.x) == smallest_u)
    {
        v = normalizeVec3(Vec3(0, -u.z, u.y, 1));
    }
    else if (abs(u.y) == smallest_u)
    {
        v = normalizeVec3(Vec3(-u.z, 0, u.x, 1));
    }
    else if (abs(u.z) == smallest_u)
    {
        v = normalizeVec3(Vec3(-u.y, u.x, 0, 1));
    }
    Vec3 w = normalizeVec3(crossProductVec3(u, v));
    double m_inverse[4][4] = {
        {u.x, v.x, w.x, 0}, {u.y, v.y, w.y, 0}, {u.z, v.z, w.z, 0}, {0, 0, 0, 1}};
    double m[4][4] = {
        {u.x, u.y, u.z, 0},
        {v.x, v.y, v.z, 0},
        {w.x, w.y, w.z, 0},
        {0, 0, 0, 1}};
    Matrix4 *M = new Matrix4(m);
    Matrix4 *M_inverse = new Matrix4(m_inverse);
    double radian = rotation->angle * M_PI / 180;
    double rotationArray[4][4] = {
        {1, 0, 0, 0},
        {0, cos(radian), -sin(radian), 0},
        {0, sin(radian), cos(radian), 0},
        {0, 0, 0, 1}};
    Matrix4 *R = new Matrix4(rotationArray);
    Matrix4 MR = multiplyMatrixWithMatrix(*M_inverse, *R);
    Matrix4 MRM = multiplyMatrixWithMatrix(MR, *M);
    matrix = multiplyMatrixWithMatrix(MRM, matrix);
}
void translateMesh(Matrix4 &matrix, const Translation *translation)
{
    double translation_[4][4] = {{1, 0, 0, translation->tx}, {0, 1, 0, translation->ty}, {0, 0, 1, translation->tz}, {0, 0, 0, 1}};
    matrix = multiplyMatrixWithMatrix(matrix, *(new Matrix4(translation_)));
}
void scaleMesh(Matrix4 &matrix, const Scaling *scaling)
{
    double scaling_[4][4] = {{scaling->sx, 0, 0, 0}, {0, scaling->sy, 0, 0}, {0, 0, scaling->sz, 0}, {0, 0, 0, 1}};
    matrix = multiplyMatrixWithMatrix(matrix, *(new Matrix4(scaling_)));
}

Matrix4 generateModelTranformation(const Mesh *mesh, const vector<Translation *> &translations, const vector<Scaling *> &scalings, const vector<Rotation *> &rotations)
{
    Matrix4 identityMatrix = getIdentityMatrix();
    for (int i = 0; i < mesh->numberOfTransformations; ++i)
    {
        char transformationType = mesh->transformationTypes[i];
        if (transformationType == 't')
        {
            translateMesh(identityMatrix, translations[mesh->transformationIds[i] - 1]);
        }
        else if (transformationType == 's')
        {
            scaleMesh(identityMatrix, scalings[mesh->transformationIds[i] - 1]);
        }
        else if (transformationType == 'r')
        {
            rotateMesh(identityMatrix, rotations[mesh->transformationIds[i] - 1]);
        }
    }
    return identityMatrix;
}
Vec3 &getCenterOfMass(const Vec3 &a, const Vec3 &b, const Vec3 &c)
{
    Vec3 *retval = new Vec3((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3, (a.z + b.z + c.z) / 3, -1);
    return *retval;
}
bool shouldCulled(const Vec4 &first, const Vec4 &second, const Vec4 &third, const Vec3 cameraPosition)
{
    Vec3 a = Vec3(first.x, second.y, third.z, -1);
    Vec3 b = Vec3(second.x, third.y, third.z, -1);
    Vec3 c = Vec3(first.x, second.y, third.z, -1);
    Vec3 normal = normalizeVec3(crossProductVec3(subtractVec3(c, b), subtractVec3(a, b)));
    Vec3 viewVector = subtractVec3(getCenterOfMass(a, b, c), cameraPosition);
    return dotProductVec3(normal, viewVector) < 0;
}
bool isVisible(double den, double num, double __TE__, double __TL__)
{
    if (den > 0)
    {
        double t = num / den;
        if (t > __TL__)
            return false;
        if (t > __TE__)
        {
            __TE__ = t;
        }
    }
    else if (den < 0)
    {
        double t = num / den;
        if (t < __TE__)
            return false;
        if (t < __TL__)
        {
            __TL__ = t;
        }
    }
    else if (num > 0)
    {
        return false;
    }
    return true;
}
bool isInView(const Camera *camera, Vec4 &v0, Vec4 &v1)
{
    double dx = v1.x - v0.x;
    double dy = v1.y - v0.y;
    double dz = v1.z - v0.z;
    double left_t = (camera->left - v0.x) / dx;
    double right_t = (camera->right - v0.x) / dx;
    double top_t = (camera->top - v0.y) / dy;
    double bottom_t = (camera->bottom - v0.y) / dy;
    double xmin = -1.0, xmax = 1.0;
    double ymin = -1.0, ymax = 1.0;
    double zmin = -1.0, zmax = 1.0;
    int __TL__ = 1, __TE__ = 0;
    bool visible = false;
    if ((dx == 0 && xmin - v0.x > 0) || (dx == 0 && v0.y - xmax > 0) || (dy == 0 && ymin - v0.y > 0) || (dy == 0 && v0.y - ymax > 0) || (dz == 0 && zmin - v0.z > 0) || (dz == 0 && v0.z - zmax > 0))
    {
        return false;
    }
    if (isVisible(dx, xmin - v0.x, __TE__, __TL__))
    {
        if (isVisible(-dx, v0.x - xmax, __TE__, __TL__))
        {
            if (isVisible(dy, ymin - v0.y, __TE__, __TL__))
            {
                if (isVisible(-dy, v0.y - ymax, __TE__, __TL__))
                {
                    if (isVisible(dz, zmin - v0.z, __TE__, __TL__))
                    {
                        if (isVisible(-dz, v0.z - zmax, __TE__, __TL__))
                        {
                            visible = true;
                            if (__TL__ < 1)
                            {
                                v1.x = v0.x + dx * __TL__;
                                v1.y = v0.y + dy * __TL__;
                                v1.z = v0.z + dz * __TL__;
                            }
                            if (__TE__ > 0)
                            {
                                v0.x = v0.x + dx * __TE__;
                                v0.y = v0.y + dy * __TE__;
                                v0.z = v0.z + dz * __TE__;
                            }
                        }
                    }
                }
            }
        }
    }
    return visible;
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
    // camera transformation
    Matrix4 *rotationMatrix = getCameraRotationMatrix(camera);
    Matrix4 *cameraTranslationMatrix = getCameraTranslationMatrix(camera);
    Matrix4 totalTransformation = multiplyMatrixWithMatrix(*rotationMatrix, *cameraTranslationMatrix);
    Matrix4 *viewPortMatrix = getViewPortMatrix(camera);
    Matrix4 *projectionMatrix = projectionTypeSwitch(camera);

    // model transformation
    for (int i = 0; i < meshes.size(); i++)
    {
        Mesh *mesh = meshes[i];
        Matrix4 modelTransformation = generateModelTranformation(mesh, translations, scalings, rotations);
        Matrix4 cameraRelativeModel = multiplyMatrixWithMatrix(*cameraTranslationMatrix, modelTransformation);
        Matrix4 projectionModel = multiplyMatrixWithMatrix(*projectionMatrix, cameraRelativeModel);
        for (int j = 0; j < mesh->triangles.size(); j++)
        {
            Triangle &triangle = mesh->triangles[j];
            Vec3 *index = vertices[triangle.getFirstVertexId() - 1];
            Vec4 first_v = Vec4(index->x, index->y, index->z, 1, index->colorId);
            index = vertices[triangle.getSecondVertexId() - 1];
            Vec4 second_v = Vec4(index->x, index->y, index->z, 1, index->colorId);
            index = vertices[triangle.getThirdVertexId() - 1];
            Vec4 third_v = Vec4(index->x, index->y, index->z, 1, index->colorId);
            first_v = multiplyMatrixWithVec4(projectionModel, first_v);
            second_v = multiplyMatrixWithVec4(projectionModel, second_v);
            third_v = multiplyMatrixWithVec4(projectionModel, third_v);
            if (cullingEnabled && shouldCulled(first_v, second_v, third_v, camera->pos))
            {
                continue;
            }
            if (isInView(camera, first_v, second_v))
            {
            }
            if (isInView(camera, second_v, third_v))
            {
            }
            if (isInView(camera, first_v, third_v))
            {
            }
        }
    }
    // clip and culling

    // rasterization process
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
    const char *str;
    XMLDocument xmlDoc;
    XMLElement *pElement;

    xmlDoc.LoadFile(xmlPath);

    XMLNode *pRoot = xmlDoc.FirstChild();

    // read background color
    pElement = pRoot->FirstChildElement("BackgroundColor");
    str = pElement->GetText();
    sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

    // read culling
    pElement = pRoot->FirstChildElement("Culling");
    if (pElement != NULL)
    {
        str = pElement->GetText();

        if (strcmp(str, "enabled") == 0)
        {
            cullingEnabled = true;
        }
        else
        {
            cullingEnabled = false;
        }
    }

    // read cameras
    pElement = pRoot->FirstChildElement("Cameras");
    XMLElement *pCamera = pElement->FirstChildElement("Camera");
    XMLElement *camElement;
    while (pCamera != NULL)
    {
        Camera *cam = new Camera();

        pCamera->QueryIntAttribute("id", &cam->cameraId);

        // read projection type
        str = pCamera->Attribute("type");

        if (strcmp(str, "orthographic") == 0)
        {
            cam->projectionType = 0;
        }
        else
        {
            cam->projectionType = 1;
        }

        camElement = pCamera->FirstChildElement("Position");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

        camElement = pCamera->FirstChildElement("Gaze");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

        camElement = pCamera->FirstChildElement("Up");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

        cam->gaze = normalizeVec3(cam->gaze);
        cam->u = crossProductVec3(cam->gaze, cam->v);
        cam->u = normalizeVec3(cam->u);

        cam->w = inverseVec3(cam->gaze);
        cam->v = crossProductVec3(cam->u, cam->gaze);
        cam->v = normalizeVec3(cam->v);

        camElement = pCamera->FirstChildElement("ImagePlane");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
               &cam->left, &cam->right, &cam->bottom, &cam->top,
               &cam->near, &cam->far, &cam->horRes, &cam->verRes);

        camElement = pCamera->FirstChildElement("OutputName");
        str = camElement->GetText();
        cam->outputFileName = string(str);

        cameras.push_back(cam);

        pCamera = pCamera->NextSiblingElement("Camera");
    }

    // read vertices
    pElement = pRoot->FirstChildElement("Vertices");
    XMLElement *pVertex = pElement->FirstChildElement("Vertex");
    int vertexId = 1;

    while (pVertex != NULL)
    {
        Vec3 *vertex = new Vec3();
        Color *color = new Color();

        vertex->colorId = vertexId;

        str = pVertex->Attribute("position");
        sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

        str = pVertex->Attribute("color");
        sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

        vertices.push_back(vertex);
        colorsOfVertices.push_back(color);

        pVertex = pVertex->NextSiblingElement("Vertex");

        vertexId++;
    }

    // read translations
    pElement = pRoot->FirstChildElement("Translations");
    XMLElement *pTranslation = pElement->FirstChildElement("Translation");
    while (pTranslation != NULL)
    {
        Translation *translation = new Translation();

        pTranslation->QueryIntAttribute("id", &translation->translationId);

        str = pTranslation->Attribute("value");
        sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

        translations.push_back(translation);

        pTranslation = pTranslation->NextSiblingElement("Translation");
    }

    // read scalings
    pElement = pRoot->FirstChildElement("Scalings");
    XMLElement *pScaling = pElement->FirstChildElement("Scaling");
    while (pScaling != NULL)
    {
        Scaling *scaling = new Scaling();

        pScaling->QueryIntAttribute("id", &scaling->scalingId);
        str = pScaling->Attribute("value");
        sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

        scalings.push_back(scaling);

        pScaling = pScaling->NextSiblingElement("Scaling");
    }

    // read rotations
    pElement = pRoot->FirstChildElement("Rotations");
    XMLElement *pRotation = pElement->FirstChildElement("Rotation");
    while (pRotation != NULL)
    {
        Rotation *rotation = new Rotation();

        pRotation->QueryIntAttribute("id", &rotation->rotationId);
        str = pRotation->Attribute("value");
        sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

        rotations.push_back(rotation);

        pRotation = pRotation->NextSiblingElement("Rotation");
    }

    // read meshes
    pElement = pRoot->FirstChildElement("Meshes");

    XMLElement *pMesh = pElement->FirstChildElement("Mesh");
    XMLElement *meshElement;
    while (pMesh != NULL)
    {
        Mesh *mesh = new Mesh();

        pMesh->QueryIntAttribute("id", &mesh->meshId);

        // read projection type
        str = pMesh->Attribute("type");

        if (strcmp(str, "wireframe") == 0)
        {
            mesh->type = 0;
        }
        else
        {
            mesh->type = 1;
        }

        // read mesh transformations
        XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
        XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

        while (pTransformation != NULL)
        {
            char transformationType;
            int transformationId;

            str = pTransformation->GetText();
            sscanf(str, "%c %d", &transformationType, &transformationId);

            mesh->transformationTypes.push_back(transformationType);
            mesh->transformationIds.push_back(transformationId);

            pTransformation = pTransformation->NextSiblingElement("Transformation");
        }

        mesh->numberOfTransformations = mesh->transformationIds.size();

        // read mesh faces
        char *row;
        char *clone_str;
        int v1, v2, v3;
        XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
        clone_str = strdup(str);

        row = strtok(clone_str, "\n");
        while (row != NULL)
        {
            int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

            if (result != EOF)
            {
                mesh->triangles.push_back(Triangle(v1, v2, v3));
            }
            row = strtok(NULL, "\n");
        }
        mesh->numberOfTriangles = mesh->triangles.size();
        meshes.push_back(mesh);

        pMesh = pMesh->NextSiblingElement("Mesh");
    }
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
    if (this->image.empty())
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            vector<Color> rowOfColors;

            for (int j = 0; j < camera->verRes; j++)
            {
                rowOfColors.push_back(this->backgroundColor);
            }

            this->image.push_back(rowOfColors);
        }
    }
    else
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            for (int j = 0; j < camera->verRes; j++)
            {
                this->image[i][j].r = this->backgroundColor.r;
                this->image[i][j].g = this->backgroundColor.g;
                this->image[i][j].b = this->backgroundColor.b;
            }
        }
    }
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
    if (value >= 255.0)
        return 255;
    if (value <= 0.0)
        return 0;
    return (int)(value);
}
void Scene::writeImageToPPMFile(Camera *camera)
{
    unsigned char *image = new unsigned char[camera->horRes * camera->verRes * 3];
    int index = 0;
    for (int j = camera->verRes - 1; j >= 0; j--)
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            image[index++] = makeBetweenZeroAnd255(this->image[i][j].r);
            image[index++] = makeBetweenZeroAnd255(this->image[i][j].g);
            image[index++] = makeBetweenZeroAnd255(this->image[i][j].b);
        }
    }
    write_ppm(camera->outputFileName.c_str(), image, camera->horRes, camera->verRes);
}

void write_ppm(const char *filename, unsigned char *data, int width, int height)
{
    FILE *outfile;

    if ((outfile = fopen(filename, "w")) == NULL)
    {
        throw std::runtime_error("Error: The ppm file cannot be opened for writing.");
    }

    (void)fprintf(outfile, "P3\n%d %d\n255\n", width, height);

    unsigned char color;
    for (size_t j = 0, idx = 0; j < height; ++j)
    {
        for (size_t i = 0; i < width; ++i)
        {
            for (size_t c = 0; c < 3; ++c, ++idx)
            {
                color = data[idx];

                if (i == width - 1 && c == 2)
                {
                    (void)fprintf(outfile, "%d", color);
                }
                else
                {
                    (void)fprintf(outfile, "%d ", color);
                }
            }
        }

        (void)fprintf(outfile, "\n");
    }

    (void)fclose(outfile);
}
/*
Original one
	Writes contents of image (Color**) into a PPM file.

void Scene::writeImageToPPMFile(Camera *camera)
{
    ofstream fout;

    fout.open(camera->outputFileName.c_str());

    fout << "P3" << endl;
    fout << "# " << camera->outputFileName << endl;
    fout << camera->horRes << " " << camera->verRes << endl;
    fout << "255" << endl;

    for (int j = camera->verRes - 1; j >= 0; j--)
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
        }
        fout << endl;
    }
    fout.close();
}
*/
/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
    string command;

    // call command on Ubuntu
    if (osType == 1)
    {
        command = "convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

    // call command on Windows
    else if (osType == 2)
    {
        command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

    // default action - don't do conversion
    else
    {
    }
}

/* 
void scale(Vec3 &v, Scaling *s)
{
    v.x *= s->sx;
    v.y *= s->sy;
    v.z *= s->sz;
}

void translate(Vec3 *v, Translation *t)
{
    v->x += t->tx;
    v->y += t->ty;
    v->z += t->tz;
}

void rotate(Vec3 *vec, Rotation *r)
{
    Vec3 *v = new Vec3(vec->x, vec->y, vec->z, vec->colorId);
    Vec3 *v2 = new Vec3(0, 0, 0, 0);
    v2->x = v->x * cos(r->angle) - v->y * sin(r->angle);
    v2->y = v->x * sin(r->angle) + v->y * cos(r->angle);
    v2->z = v->z;
    vec->x = v2->x;
    vec->y = v2->y;
    vec->z = v2->z;
}

void modelingTransformation(vector<Vec3 *> vertices, vector<Scaling *> scalings, vector<Rotation *> rotations, vector<Translation *> translations)
{
    // TODO : Implement this function
} */