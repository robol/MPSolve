#include <mps/mps.h>
#include <gtk/gtk.h>

G_BEGIN_DECLS

#if GTK_MAJOR_VERSION < 3

#ifndef gtk_widget_get_allocated_width
#define gtk_widget_get_allocated_width(widget) (widget->allocation.width)
#endif

#ifndef gtk_widget_get_allocated_height
#define gtk_widget_get_allocated_height(widget) (widget->allocation.height)
#endif

#endif

#ifndef MPS_ITERATION_LOGGER_H_
#define MPS_ITERATION_LOGGER_H_

#define MPS_TYPE_ITERATION_LOGGER               (mps_iteration_logger_get_type ())
#define MPS_ITERATION_LOGGER(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), MPS_TYPE_ITERATION_LOGGER, MpsIterationLogger))
#define MPS_IS_ITERATION_LOGGER(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), MPS_TYPE_ITERATION_LOGGER))
#define MPS_ITERATION_LOGGER_CLASS(klass)       (G_TYPE_CHECK_CLASS_CAST ((klass), MPS_TYPE_ITERATION_LOGGER, MpsIterationLoggerClass))
#define MAMAN_IS_BAR_CLASS(klass)               (G_TYPE_CHECK_CLASS_TYPE ((klass), MPS_TYPE_ITERATION_LOGGER))
#define MPS_ITERATION_LOGGER_GET_CLASS(obj)     (G_TYPE_INSTANCE_GET_CLASS ((obj), MPS_TYPE_ITERATION_LOGGER, MpsIterationLoggerClass)))

typedef struct _MpsIterationLogger              MpsIterationLogger;
typedef struct _MpsIterationLoggerClass         MpsIterationLoggerClass;

struct _MpsIterationLogger {
        GtkWindow parent_instance;

        /* <private declarations>*/
        GtkWidget * drawing_area;
        guint timeout_source;
        mps_context * ctx;

        gboolean drawing;
        
        /* Scale of the plot */
        double x_scale;
        double y_scale;

        mps_approximation ** approximations;
        int degree;

        /* Handling of the zomming process */
        gboolean zooming;
        gint     zoom_rect_x;
        gint     zoom_rect_y;

        double real_center;
        double imag_center;

        gboolean exit;

        pthread_mutex_t *drawing_lock;
};

struct _MpsIterationLoggerClass {
        GtkWindowClass parent_class;
};

GType mps_iteration_logger_get_type (void);

/**
 * @brief Allocate a new MpsIterationLogger. 
 */
MpsIterationLogger* mps_iteration_logger_new (void);

/**
 * @brief Assing an mps_context to this MpsIterationLogger, so it will display the
 * approximations as soon as they change. 
 */
void mps_iteration_logger_set_mps_context (MpsIterationLogger * logger, mps_context * context);

/**
 * @brief Assign to the iteration logger a set of static roots to be displayed.
 */
void mps_iteration_logger_set_roots (MpsIterationLogger * logger, mps_approximation ** approximations, int degree);


#endif

G_END_DECLS


