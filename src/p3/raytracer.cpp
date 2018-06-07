/**
 * @file raytacer.cpp
 * @brief Raytracer class
 *
 * 
 *
 * @author HingOn Miu (hmiu)
 * 
 */

#include "raytracer.hpp"
#include "scene/scene.hpp"

#include <SDL_timer.h>
#include <iostream>
#include <random>

#ifdef OPENMP // just a defense in case OpenMP is not installed.

#include <omp.h>

#endif
namespace _462 {

// max number of threads OpenMP can use. Change this if you like.
#define MAX_THREADS 8
// since the average sunlight in the air is yellow-orange color
// so i choose its wavelength as the wavelength of the eye ray
// shooting from camera in the air
#define SUNLIGHT_WAVELENGTH 598
// the number of shadow ray to sample the spherical light
#define NUM_SHADOWRAY 5

static const unsigned STEP_SIZE = 8;
// the vertical and horizontal vectors for the camera len plane
Vector3 len_vertical_scale;
Vector3 len_horizontal_scale;
// distance between the focus point and the screen
real_t focal_length;

Raytracer::Raytracer()
    : scene(0), width(0), height(0) { }

// random real_t in [0, 1)
static inline real_t random()
{
    return real_t(rand())/RAND_MAX;
}

Raytracer::~Raytracer() { }

/**
 * Initializes the raytracer for the given scene. Overrides any previous
 * initializations. May be invoked before a previous raytrace completes.
 * @param scene The scene to raytrace.
 * @param width The width of the image being raytraced.
 * @param height The height of the image being raytraced.
 * @return true on success, false on error. The raytrace will abort if
 *  false is returned.
 */
bool Raytracer::initialize(Scene* scene, size_t num_samples,
			   size_t width, size_t height)
{
    /*
     * omp_set_num_threads sets the maximum number of threads OpenMP will
     * use at once.
     */
#ifdef OPENMP
    omp_set_num_threads(MAX_THREADS);
#endif
    this->scene = scene;
    this->num_samples = num_samples;
    this->width = width;
    this->height = height;

    current_row = 0;

    Ray::init(scene->camera);
    scene->initialize();

    const SphereLight* lights = scene->get_lights();

    // TODO any initialization or precompuation before the trace
    
    // get the distance between sceen and focus point
    focal_length = 0.1 * scene->camera.get_far_clip();
    
    return true;
}

/**
 * Set up the camera len plane's vertical and horizontal scale.
 * for depth of field effect with updated camera info
 * @param scene The scene to raytrace.
 * @return void
 */
void setup_depth_of_field(const Scene* scene) {
    if (scene->camera_len_length != 0) {
        Vector3 c = scene->camera.get_position();
        Vector3 n = scene->camera.get_direction();
        Vector3 up = scene->camera.get_up();
        // the length is measured here in 0.1 per unit
        real_t length = ((real_t)scene->camera_len_length)/10.0;
        // determine the vertical and horizontal vectors for the len plane
        len_vertical_scale = length * up;
        len_horizontal_scale = length * normalize(cross(up, n));
    }
    return;
}

/**
 * Cast a ray to the scene and check if any object is hit.
 * @param scene The scene to raytrace.
 * @param r The ray to cast.
 * @param inter The details of intersected object.
 * @return true if hit an object, and false otherwise.
 */
bool cast_ray(const Scene* scene, const Ray& r, Intersection& inter) {
    size_t i;
    // record the closest t found
    real_t closest_found_t = t1;
    // check if the ray does hit an object 
    for (i = 0; i < scene->num_geometries(); i++) {
        // instead of comparing distance here inefficiently,
        // intersection_test check closes_found_t before editing
        // inter, and so inter will be the closest intersection
        // after the for loop.
        ((scene->get_geometries())[i])->
            intersection_test(r, inter, &closest_found_t);
    }
    
    // check if the closest_found_t is still the upper time bound
    // to confirm whether there is an intersection
    return closest_found_t > t0 && closest_found_t < t1;
}

/**
 * Use Monte Carlo method to sample the spherical light.
 * @param light The spherical light.
 * @return The new light position.
 */
Vector3 monte_carlo_sampling(const SphereLight light) {
    real_t x = random_gaussian();
    real_t y = random_gaussian();
    real_t z = random_gaussian();
    // generate a random point on a unit sphere
    Vector3 surface_p = normalize(Vector3(x,y,z));
    // scale the point by sphere radius and transform it by sphere center 
    return surface_p * light.radius + light.position;   
}

/**
 * Use Blinn-Phong illumination to calculate direct illumination.
 * @param scene The scene to render.
 * @param r The ray that hits an object.
 * @param inter The details of the intersected object.
 * @return The color of the intersection.
 */
Color3 direct_illumination(const Scene *scene, const Ray &r, 
                           Intersection& inter) {
    // support the basic Blinn-Phong illumination
    // compute the ambient term
    Color3 a = scene->ambient_light * inter.material.ambient;
    
    size_t i, j;
    real_t dis, non_blocked;
    Vector3 dir, volume_light_dir;
    Color3 d = Color3::Black();
    Color3 light_color;
    real_t N_L;
    Intersection tmp;
    Ray light_r;
    // compute the diffuse term
    for (i = 0; i < scene->num_lights(); i++) {
        // light direction
        dir = scene->get_lights()[i].position - 
              inter.world_intersection;
        // distance from the light
        dis = length(dir);
	
	// record the number of non-blocked shadow rays
        non_blocked = NUM_SHADOWRAY;
        for (j = 0; j < NUM_SHADOWRAY; j++) {     
            // use monte carlo method to sample volumne light source position
            volume_light_dir = monte_carlo_sampling(scene->get_lights()[i])
                               - inter.world_intersection; 
            light_r = Ray(inter.world_intersection, volume_light_dir); 
            // check the shadow, whether the light is blocked by closer object
            if (cast_ray(scene, light_r, tmp) == true &&
                dis > length(tmp.world_intersection-inter.world_intersection)){
                non_blocked--;
            }
        } 
        
        // compute the color of the light with attenuation
        light_color = (scene->get_lights()[i]).color *  
             (1.0 / (scene->get_lights()[i].attenuation.constant +
             scene->get_lights()[i].attenuation.linear * dis +
             scene->get_lights()[i].attenuation.quadratic * dis * dis));       
        // dot product of intersection normal and light vector
        N_L = dot(inter.world_normal, normalize(dir));
                  //inter.world_intersection - 
        // max(N_L, 0.0)
        if (N_L < 0.0) {
            N_L = 0.0;
        }        
        
        d += (non_blocked/(real_t)NUM_SHADOWRAY) * light_color * 
             inter.material.diffuse * N_L;
    }
    
    return inter.texture_color * (a + d);
}

/**
 * Calculate the reflected ray.
 * @param e The intersection location.
 * @param d The direction of incoming ray.
 * @param n The normal vector on surface.
 * @return The refected ray.
 */
Ray reflect(Vector3 e, Vector3 d, Vector3 n) {
    // reflected_dir is d - 2*(d.n)*n
    Vector3 reflect_dir = d - 2.0 * dot(d, n) * n;
    return Ray(e, reflect_dir);
}

/**
 * Check whether there is refracted ray from the intersected object.
 * @param d The direction of incoming ray.
 * @param norm The normal vector of surface.
 * @param n The current medium's refractive index.
 * @param n The intersected medium's refractive index.
 * @param t The refracted ray direction.
 * @return true if there is refracted ray, and false otherwise.
 */
bool refract(Vector3 d, Vector3 norm, real_t n, real_t nt, Vector3* t) {
    // if the number under the square root is negative, there is no
    // refracted ray and all of the energy is reflected, which is
    // known as total internal reflection.
    real_t sq_rt = 1.0 - (n * n * (1.0 - 
                   (dot(d, norm) * dot(d, norm))))/ (nt * nt);

    if (sq_rt < 0.0) {
        return false;
    }
    else {
        Vector3 tmp = 
            (n * (d - norm * dot(d, norm))) / nt - norm * sqrt(sq_rt);
        t->x = tmp.x;
        t->y = tmp.y;
        t->z = tmp.z;
        return true;
    }
}

/**
 * Get the current and next refractive index and wavelength from stack.
 * @param scene The scene to render.
 * @param inter The details of intersection.
 * @param entering Indicates whether entering or leaving dielectric.
 * @param w The current medium's traveling ray's wavelength.
 * @param wt The next medium's traveling ray's wavelength.
 * @param n The current medium's refractive index.
 * @param nt The next medium's refractive index.
 * @param refractive_index_stack The stack to hold all refractive index.
 * @param wavelength_stack The stack to hold all wavelengths.
 * @param i The current position on the stack.
 * @return void.
 */
void access_stack(const Scene *scene, Intersection& inter,
                  bool entering, real_t *w, real_t* wt, real_t *n, real_t *nt,
                  std::vector<real_t>& refractive_index_stack,
                  std::vector<real_t>& wavelength_stack, int i) {
    if (entering) {
        *n = refractive_index_stack[i];
        *nt = inter.material.refractive_index;
        // add the intersected object's refractive_index to stack
        refractive_index_stack.push_back(*nt);
        // compute the intersected object's traveling light wavelength
        if (scene->Beers_Law) {
            // according to Snell's Law, two medium's infractive index ratio
            // inversely equals two medium's wavelength ratio, such that
            // n1/n2 = w2/w1, and so new wavelength = (n1/n2) * w1
            *w = wavelength_stack[i];
            *wt = ((*n)/(*nt)) * (*w);
            // add the intersected object's wavelength to stack
            wavelength_stack.push_back(*wt);
        }
    }
    else {
        *n = inter.material.refractive_index;
        int j = (i > 0) ? (i - 1) : 0;
        *nt = refractive_index_stack[j];
        if (scene->Beers_Law) {
            *wt = wavelength_stack[j];
            *w = (*wt) / ((*n)/(*nt)); 
        }
    }
}

/**
 * Calculate attenuation factor for Beer's Law equation.
 * @param t The interection distance.
 * @param n The current medium's refractive index.
 * @param w The current medium's traveling ray's wavelength.
 * @param The attenuation factor.
 */
real_t get_attenuation_factor(real_t t, real_t n, real_t w) { 
    // attenuation factor k is calculated as exp(-a * t) in shirley textbook
    // where a is the attenuation coefficient and t is the traveled distance.
    // according to Beer-Lambert Law on wikipedia, the calculation of the
    // attenuation coefficient a is (4 * pi * n)/w where n is the refractive
    // index and w is the traveling light's wavelength.
    real_t a = (4.0 * PI * n)/w;
    real_t k = exp(-a * t);
    return k;
}

/**
 * Get the pixel color by calculating direct illumination, reflection,
 * and refraction.
 * @param scene The scene to render.
 * @param r The ray to cast.
 * @param depth The recursion depth.
 * @param refractive_index_stack The stack to hold all refractive index.
 * @param wavelength_stack The stack to hold all wavelengths.
 * @param current_stack_position The current position on the stack.
 * @return The color of the pixel.
 */
Color3 get_pixel_color(const Scene *scene, const Ray& r, 
                       int depth, std::vector<real_t>& refractive_index_stack,
                       std::vector<real_t>& wavelength_stack,
                       int current_stack_position) {
    // the maximum recursive depth for tracing the ray from screen
    if (depth > MAX_DEPTH) {
         return Color3::Black();
    }
    
    Intersection inter; 
    // nothing is hit by the eye ray 
    if (cast_ray(scene, r, inter) == false) {
        return scene->background_color;
    }
    
    // compute the reflected ray 
    Ray reflect_r = reflect(inter.world_intersection, normalize(r.d), 
                            inter.world_normal); 

    // handle opaque objects
    if (inter.material.refractive_index == 0.0) {
        // compute the color by direct illumination
        Color3 direct_color = direct_illumination(scene, r, inter);
        // sum the color from direct illumination and 
        // color from specular term's reflected ray computation
        return direct_color + inter.texture_color * inter.material.specular *
               get_pixel_color(scene, reflect_r, depth + 1,
                               refractive_index_stack, 
                               wavelength_stack,
                               current_stack_position);
    }
    // handle dielectrics
    else {
        real_t c, k, w, wt, n, nt;
        Vector3 t;
        bool entering;
        
        // entering the dielectric 
        if (dot(r.d, inter.world_normal) < 0.0) {
            entering = true;
            access_stack(scene, inter, entering, &w, &wt, &n, &nt,
                         refractive_index_stack, wavelength_stack,
                         current_stack_position);
            k = 1.0;
            refract(normalize(r.d), inter.world_normal, n, nt, &t);
            c = dot(-normalize(r.d), inter.world_normal);

        }
        // leaving the dialectric
        else {
            entering = false;
            access_stack(scene, inter, entering, &w, &wt, &n, &nt,
                         refractive_index_stack, wavelength_stack,
                         current_stack_position);
            // attenuation factor according to Beer's Law
            k = (scene->Beers_Law) ? 
                get_attenuation_factor(inter.local_t, n, w) : 1.0;
            if (refract(normalize(r.d), -inter.world_normal, n, nt, &t)) {
                c = dot(normalize(t), inter.world_normal); 
            }
            // total internal reflection
            else {
                // R = 1 for total internal reflection
                return k * get_pixel_color(scene, reflect_r, depth + 1,
                                           refractive_index_stack,
                                           wavelength_stack,
                                           current_stack_position);
            }
        }
        
        // compute Fresnel coefficient R
        real_t R0 = (real_t)pow((double)((nt - 1.0)/(nt + 1.0)), 2.0);
        real_t R = R0 + (1.0 - R0) * (real_t)pow(1.0 - c, 5.0);
        Ray refract_r = Ray(inter.world_intersection, normalize(t));
        
        // schlick approximation of the Fresnel effect
        // point color = k * (R * reflection color + (1-R) * refraction color)
        Color3 reflection_term = 
            R * get_pixel_color(scene, reflect_r, depth + 1, 
            refractive_index_stack, wavelength_stack, current_stack_position);
        
        // only edit the stack position  when entering or leaving an 
        // object happens, and so reflecting outside or inside an object 
        // would not change the current stack position
        if (entering) {
            // entering object means the new current refractive index
            // would be the one next on stack
            return k * (reflection_term +
                (1.0 - R) * get_pixel_color(scene, refract_r, depth + 1,
                refractive_index_stack, wavelength_stack,
                current_stack_position + 1));
        }
        else {
            // leaving object means the new current refractive index
            // would be the one previous on stack
            return k * (reflection_term +
                (1.0 - R) * get_pixel_color(scene, refract_r, depth + 1,
                refractive_index_stack, wavelength_stack,
                current_stack_position - 1));
        }

    }
}

/**
 * Get the camera eye position.
 * @param scene The scene to render.
 * @return The position of camera eye.
 */
Vector3 get_camera_eye(const Scene* scene) {
    // user does not specify the camera len length, so the
    // depth of field effect is disable and simply return
    // the fixed camera position
    if (scene->camera_len_length == 0) {
        return scene->camera.get_position();
    }

    Vector3 e = scene->camera.get_position();
    // shift the camera eye vertically and horizontally with the len
    e += (len_vertical_scale - random() * len_vertical_scale * 2.0);
    e += (len_horizontal_scale - random() * len_horizontal_scale * 2.0);
    return e;
}

/**
 * Get the focus point for the ray from the len.
 * @param scene The scene to render.
 * @param x The x-coordinate of the pixel to trace.
 * @param y The y-coordinate of the pixel to trace.
 * @param width The width of the screen in pixels.
 * @param height The height of the screen in pixels.
 * @return The location of the focus point.
 */
Vector3 get_focus_point(const Scene*scene, size_t x, size_t y,
			size_t width, size_t height) {
    // concept reference for depth of field effect is from
    // http://http.developer.nvidia.com/GPUGems/gpugems_ch23.html
    real_t dx = real_t(1)/width;
    real_t dy = real_t(1)/height;

    // get the direction of ray shooting from center of the pixel
    real_t center_i = real_t(2)*(real_t(x) + 0.5)*dx - real_t(1);
    real_t center_j = real_t(2)*(real_t(y) + 0.5)*dy - real_t(1);
    Vector3 focus_dir = Ray::get_pixel_dir(center_i, center_j);
    // compute the location of the focus point
    return scene->camera.get_position() + focal_length * focus_dir;
}

/**
 * Performs a raytrace on the given pixel on the current scene.
 * The pixel is relative to the bottom-left corner of the image.
 * @param scene The scene to trace.
 * @param x The x-coordinate of the pixel to trace.
 * @param y The y-coordinate of the pixel to trace.
 * @param width The width of the screen in pixels.
 * @param height The height of the screen in pixels.
 * @return The color of that pixel in the final image.
 */
Color3 Raytracer::trace_pixel(const Scene* scene,
			      size_t x,
			      size_t y,
			      size_t width,
			      size_t height)
{
    assert(x < width);
    assert(y < height);
    
    real_t dx = real_t(1)/width;
    real_t dy = real_t(1)/height;   
    // compute the location of the focus point
    Vector3 focus = get_focus_point(scene, x, y, width, height);
    Color3 res = Color3::Black();
    unsigned int iter;
    for (iter = 0; iter < num_samples; iter++)
    {
        // compute the ray shooting our from the pixel
        Vector3 e = get_camera_eye(scene); 
        Vector3 dir;
        if (scene->camera_len_length != 0){
            // depth of field effect is turned on
            dir = normalize(focus - e);
        }
        else {
            // depth of field effect is turned off
            real_t i = real_t(2)*(real_t(x) + random())*dx - real_t(1);
            real_t j = real_t(2)*(real_t(y) + random())*dy - real_t(1);
            dir = Ray::get_pixel_dir(i, j);
        }
        Ray r = Ray(e, dir);
        
        // record all the refractive index to calculate the Frensel equation
        std::vector<real_t> refractive_index_stack;
        // place the air's refractive index to stack to start with
        refractive_index_stack.push_back(scene->refractive_index); 
        // record all the wavelengths to calculate the Beer's Law
        std::vector<real_t> wavelength_stack;
        if (scene->Beers_Law) {
            // place the air's sunlight wavelength to stack to start with
            wavelength_stack.push_back((real_t)SUNLIGHT_WAVELENGTH); 
        }
        res += get_pixel_color(scene, r, 0, refractive_index_stack, 
                               wavelength_stack, 0);

    }
    return res*(real_t(1)/num_samples);
}

/**
 * Raytraces some portion of the scene. Should raytrace for about
 * max_time duration and then return, even if the raytrace is not copmlete.
 * The results should be placed in the given buffer.
 * @param buffer The buffer into which to place the color data. It is
 *  32-bit RGBA (4 bytes per pixel), in row-major order.
 * @param max_time, If non-null, the maximum suggested time this
 *  function raytrace before returning, in seconds. If null, the raytrace
 *  should run to completion.
 * @return true if the raytrace is complete, false if there is more
 *  work to be done.
 */
bool Raytracer::raytrace(unsigned char* buffer, real_t* max_time)
{
    // TODO Add any modifications to this algorithm, if needed.
    
    // update the len plane orientation and location with updated camera info
    setup_depth_of_field(scene);
    
    static const size_t PRINT_INTERVAL = 64;

    // the time in milliseconds that we should stop
    unsigned int end_time = 0;
    bool is_done;

    if (max_time)
    {
        // convert duration to milliseconds
        unsigned int duration = (unsigned int) (*max_time * 1000);
        end_time = SDL_GetTicks() + duration;
    }

    // until time is up, run the raytrace. we render an entire group of
    // rows at once for simplicity and efficiency.
    for (; !max_time || end_time > SDL_GetTicks(); current_row += STEP_SIZE)
    {
        // we're done if we finish the last row
        is_done = current_row >= height;
        // break if we finish
        if (is_done) break;

        int loop_upper = std::min(current_row + STEP_SIZE, height);

        // This tells OpenMP that this loop can be parallelized.
#pragma omp parallel for
        for (int c_row = current_row; c_row < loop_upper; c_row++)
        {
            /*
             * This defines a critical region of code that should be
             * executed sequentially.
             */
#pragma omp critical
            {
                if (c_row % PRINT_INTERVAL == 0)
                    printf("Raytracing (Row %d)\n", c_row);
            }

            for (size_t x = 0; x < width; x++)
            {
                // trace a pixel
                Color3 color = trace_pixel(scene, x, c_row, width, height);
                // write the result to the buffer, always use 1.0 as the alpha
                color.to_array(&buffer[4 * (c_row * width + x)]);
            }
        }
    }

    if (is_done) printf("Done raytracing!\n");

    return is_done;
}

} /* _462 */
