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
        GObject parent_instance;

        /* <private declarations>*/
        GtkWidget * window;
        GtkWidget * drawing_area;
        mps_context * ctx;
        
        double x_scale;
        double y_scale;

        gboolean exit;
};

struct _MpsIterationLoggerClass {
        GObjectClass parent_class;
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


#endif

G_END_DECLS


