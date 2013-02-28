#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_GTK

#define _MPS_PRIVATE
#include <iteration-logger.h>
#include <mps/mps.h>
#include <math.h>

G_DEFINE_TYPE (MpsIterationLogger, mps_iteration_logger, G_TYPE_OBJECT);

static void mps_iteration_logger_build_interface (MpsIterationLogger*);

static void
mps_iteration_logger_dispose (GObject *object)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (object);

  logger->exit = TRUE;
  logger->ctx = NULL;
}

static void
mps_iteration_logger_finalize (GObject *object)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (object);
  g_object_unref (logger->window);
}

static void
mps_iteration_logger_class_init (MpsIterationLoggerClass *klass)
{
  GObjectClass * gobject_class = G_OBJECT_CLASS (klass);

  gobject_class->dispose = mps_iteration_logger_dispose;
  gobject_class->finalize = mps_iteration_logger_finalize;
}

static void
mps_iteration_logger_init (MpsIterationLogger * logger)
{
  logger->ctx = NULL;
  logger->window = NULL;
  logger->drawing_area = NULL;
  logger->exit = FALSE;

  /* Setup neutral zomming */
  logger->zooming = false;
  logger->real_center = 0.0;
  logger->imag_center = 0.0;
  logger->x_scale = logger->y_scale = 1.0;

  mps_iteration_logger_build_interface (logger);
}

MpsIterationLogger*
mps_iteration_logger_new ()
{
  MpsIterationLogger* logger = MPS_ITERATION_LOGGER (g_object_new (MPS_TYPE_ITERATION_LOGGER, 0, NULL));
  return logger;
}

void
mps_iteration_logger_set_mps_context (MpsIterationLogger * logger, mps_context * context)
{
  logger->ctx = context;
}

static void mps_iteration_logger_on_drawing_area_draw (GtkWidget* , cairo_t*, MpsIterationLogger * logger);

#if GTK_MAJOR_VERSION < 3
static void
mps_iteration_logger_on_expose_event (GtkWidget * widget, GdkEvent * event, MpsIterationLogger * logger)
{
  cairo_t * cr = gdk_cairo_create (widget->window);
  mps_iteration_logger_on_drawing_area_draw (widget, cr, logger);
}
#endif

#if GTK_MAJOR_VERSION >= 3
static void
get_pointer_position (GtkWidget *widget,
                      gint *x,
                     gint *y)
{
 GdkWindow *window;
 GdkDisplay *display;
 GdkDeviceManager *device_manager;
 GdkDevice *device;

 window = gtk_widget_get_window (widget);
 display = gdk_window_get_display (window);
 device_manager = gdk_display_get_device_manager (display);
 device = gdk_device_manager_get_client_pointer (device_manager);

 gdk_window_get_device_position (window, device, x, y, NULL);
}
#else
#define get_pointer_position(text_view,x,y) gtk_widget_get_pointer(text_view, x, y)
#endif

static void
mps_iteration_logger_on_drawing_area_draw (GtkWidget * widget,
                                           cairo_t * cr, MpsIterationLogger * logger)
{
  int width, height;
  mps_operation operation = logger->ctx->operation;

  if (!logger->ctx || logger->exit)
    return;

  int degree = mps_context_get_degree (logger->ctx);

  width = gtk_widget_get_allocated_width (widget);
  height = gtk_widget_get_allocated_height (widget);

  /* Draw background */
  cairo_rectangle (cr, 0, 0, width, height);
  cairo_set_source_rgb (cr, 1, 1, 1);
  cairo_fill (cr);

  /* Draw x and y axes */
  int x_level = (logger->imag_center / logger->y_scale + 1) * height * 0.5;
  int y_level = (-logger->real_center / logger->x_scale + 1) * width * 0.5;

  cairo_set_line_width (cr, 1);
  cairo_set_source_rgb (cr, 0, 0, 0);

  cairo_move_to (cr, 6, x_level);
  cairo_line_to (cr, width - 6, x_level);

  cairo_move_to (cr, y_level, height - 6);
  cairo_line_to (cr, y_level, 6);

  cairo_stroke (cr);

  /* Draw points if present */
  if (logger->ctx->root && (operation == MPS_OPERATION_REGENERATION ||
      operation == MPS_OPERATION_ABERTH_FP_ITERATIONS ||
      operation == MPS_OPERATION_ABERTH_DPE_ITERATIONS ||
      operation == MPS_OPERATION_ABERTH_MP_ITERATIONS))
    {
      int i;
      double x, y;
      cairo_set_source_rgb (cr, 0.9, 0.1, 0.1);

#define PADDING (24)
      /* Check if we can draw in here. */
      if (width < 2 * PADDING || height < 2 * PADDING)
        return;

      for (i = 0; i < degree; i++)
        {
          x =  (cplx_Re (logger->ctx->root[i]->fvalue) - logger->real_center) *
               (0.5 / logger->x_scale * width - PADDING) + width / 2;
          y =  -(cplx_Im (logger->ctx->root[i]->fvalue) - logger->imag_center) * 
               (0.5 / logger->y_scale * height - PADDING) + height / 2;
          cairo_arc (cr, x, y, 1.3, 0, 6.29);
          cairo_fill (cr);
        }
#undef PADDING

      if (logger->zooming)
        {
          int x, y;
          get_pointer_position (GTK_WIDGET (logger->drawing_area), &x, &y);

          cairo_set_source_rgb (cr, 0.4, 0.7, 0.7);
          cairo_rectangle (cr, 
              logger->zoom_rect_x, logger->zoom_rect_y, 
              x - logger->zoom_rect_x, y - logger->zoom_rect_y);
          cairo_stroke (cr);

          cairo_set_source_rgba (cr, 0.7, 1.0, 0.7, 0.3);
          cairo_rectangle (cr, 
              logger->zoom_rect_x, logger->zoom_rect_y, 
              x - logger->zoom_rect_x, y - logger->zoom_rect_y);
          cairo_fill (cr);
        }

    }
}

static void
mps_iteration_logger_on_delete_event (GtkWidget * widget, GdkEvent* event, void * user_data)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (user_data);
  logger->exit = true;
}

static gboolean 
mps_iteration_logger_update_drawing_area (gpointer user_data)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (user_data);

  if (!logger->exit)
    gtk_widget_queue_draw (logger->drawing_area);
  else
  {
    gtk_widget_destroy (logger->window);
  }

  return !logger->exit;
}

static void
mps_iteration_logger_on_zoom_out_clicked (GtkWidget * button, gpointer user_data)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (user_data);

  logger->x_scale *= 2;
  logger->y_scale *= 2;
}

static void
mps_iteration_logger_on_da_button_press (GtkWidget *widget, GdkEventButton * event, gpointer user_data)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (user_data);

  int width = gtk_widget_get_allocated_width (widget);
  int height = gtk_widget_get_allocated_height (widget);

#define HOR_COORDS_TO_POINTS(x) ((x * 2.0 / width) - 1.0)
#define VER_COORDS_TO_POINTS(y) ((-y * 2.0 / height) + 1.0)

  if (!logger->zooming)
    {
      logger->zooming = true;
      logger->zoom_rect_x = event->x;
      logger->zoom_rect_y = event->y;
    }
  else 
    {
      logger->zooming = FALSE;

      logger->real_center = HOR_COORDS_TO_POINTS ((event->x + logger->zoom_rect_x) / 2);
      logger->imag_center = VER_COORDS_TO_POINTS ((event->y + logger->zoom_rect_y) / 2);

      logger->x_scale *= (fabs (1.0 * event->x - logger->zoom_rect_x) / width);
      logger->y_scale *= (fabs (1.0 * event->y - logger->zoom_rect_y) / height);
    }
}

static void
mps_iteration_logger_reset_zoom (MpsIterationLogger * logger)
{
  logger->zooming = FALSE;
  logger->real_center = 0.0;
  logger->imag_center = 0.0;
  logger->x_scale = logger->y_scale = 1.0;
}

static void
mps_iteration_logger_on_reset_zoom_clicked (GtkWidget * button, gpointer user_data)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (user_data);
  mps_iteration_logger_reset_zoom (logger);
}

static void
mps_iteration_logger_build_interface (MpsIterationLogger * logger)
{
  /* The window we are using for this simple example, where all the widgets will
   * be packed. */
  logger->window = gtk_window_new (GTK_WINDOW_TOPLEVEL);

  /* We will use this drawing area to draw the computed roots. */
  logger->drawing_area = gtk_drawing_area_new ();

  /* Setup the drawing area to to handle the appropriate events */
  gtk_widget_add_events (logger->drawing_area, GDK_BUTTON_PRESS_MASK);
  g_signal_connect (logger->drawing_area, "button-press-event", 
    G_CALLBACK (mps_iteration_logger_on_da_button_press), logger);

  GtkWidget *box = gtk_box_new (GTK_ORIENTATION_VERTICAL, 6);

  gtk_box_pack_start (GTK_BOX (box), logger->drawing_area, TRUE, TRUE, 0);

  /* Create the range adjuster */
  GtkWidget * zoom_out = gtk_button_new_with_label ("Zoom out");
  GtkWidget * zoom_in  = gtk_button_new_with_label ("Reset zoom");

  g_signal_connect (zoom_out, "clicked", G_CALLBACK (mps_iteration_logger_on_zoom_out_clicked),
    logger);
  g_signal_connect (zoom_in, "clicked", G_CALLBACK (mps_iteration_logger_on_reset_zoom_clicked),
    logger);

  GtkWidget * zoom_box = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 6);

  gtk_box_pack_start (GTK_BOX (zoom_box), zoom_in, TRUE, TRUE, 0);
  gtk_box_pack_start (GTK_BOX (zoom_box), zoom_out, TRUE, TRUE, 0);

  gtk_box_pack_start (GTK_BOX (box), zoom_box, FALSE, TRUE, 0);

  gtk_container_add (GTK_CONTAINER (logger->window), box);

#if GTK_MAJOR_VERSION < 3
  g_signal_connect (logger->drawing_area, "expose-event", 
    G_CALLBACK (mps_iteration_logger_on_expose_event), logger);
#else
  g_signal_connect (logger->drawing_area, "draw", 
    G_CALLBACK (mps_iteration_logger_on_drawing_area_draw), logger);
#endif

  g_signal_connect (logger->window, "delete_event", 
    G_CALLBACK (mps_iteration_logger_on_delete_event), logger);

  /* Start the loop */
  g_timeout_add (1000.0 / 25, mps_iteration_logger_update_drawing_area, logger);

  gtk_widget_set_size_request (logger->window, 640, 640);
}

#endif