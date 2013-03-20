#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_GRAPHICAL_DEBUGGER

#define _MPS_PRIVATE
#include <iteration-logger.h>
#include <mps/mps.h>
#include <math.h>

G_DEFINE_TYPE (MpsIterationLogger, mps_iteration_logger, GTK_TYPE_WINDOW);

static void mps_iteration_logger_build_interface (MpsIterationLogger*);

static void
mps_iteration_logger_dispose (GObject *object)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (object);
  if (logger->timeout_source <= 0 || !g_source_remove (logger->timeout_source))
    g_warning ("Source not found");
  logger->exit = TRUE;
}

static void
mps_iteration_logger_finalize (GObject *object)
{
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
  logger->drawing_area = NULL;
  logger->exit = FALSE;

  /* Setup neutral zomming */
  logger->zooming = false;
  logger->real_center = 0.0;
  logger->imag_center = 0.0;
  logger->x_scale = logger->y_scale = 1.0;

  logger->degree = 0;
  logger->approximations = NULL;

  logger->drawing = FALSE;

  mps_iteration_logger_build_interface (logger);
}

MpsIterationLogger*
mps_iteration_logger_new ()
{
  MpsIterationLogger* logger = MPS_ITERATION_LOGGER (g_object_new (MPS_TYPE_ITERATION_LOGGER, "type",
    GTK_WINDOW_TOPLEVEL, NULL));
  return logger;
}

void
mps_iteration_logger_set_mps_context (MpsIterationLogger * logger, mps_context * context)
{
  logger->ctx = context;
  // mps_iteration_logger_set_roots (logger, logger->ctx->root, mps_context_get_degree (context));
}

void
mps_iteration_logger_set_roots (MpsIterationLogger * logger, mps_approximation ** approximations, int n)
{
  logger->approximations = approximations;
  logger->degree = n;
  logger->ctx = NULL;
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

static gdouble
mps_iteration_logger_x_coords_to_points (MpsIterationLogger * logger, gdouble x)
{
  int width = gtk_widget_get_allocated_width (logger->drawing_area);
  return ((x * 2.0 / width) - 1.0) * logger->x_scale + logger->real_center;
}

static gdouble
mps_iteration_logger_y_coords_to_points (MpsIterationLogger * logger, gdouble y)
{
  int height = gtk_widget_get_allocated_height (logger->drawing_area);
  return 2.0 * (-y + height / 2) / height * logger->y_scale + logger->imag_center;
}

static gdouble
mps_iteration_logger_x_points_to_coords (MpsIterationLogger * logger, gdouble x)
{
  int width = gtk_widget_get_allocated_width (logger->drawing_area);
  return (x - logger->real_center) * (0.5 / logger->x_scale * width) + width / 2;
}

static gdouble
mps_iteration_logger_y_points_to_coords (MpsIterationLogger * logger, gdouble y)
{
  int height = gtk_widget_get_allocated_height (logger->drawing_area);
  return -(y - logger->imag_center) * (0.5 / logger->y_scale * height) + height / 2;
}

static void
mps_iteration_logger_draw_x_ticks (MpsIterationLogger * logger, cairo_t * cr)
{
  int height = gtk_widget_get_allocated_height (logger->drawing_area);
  cairo_text_extents_t extents;

  cairo_set_source_rgb (cr, 0.4, 0.4, 0.4);

  /* Get a proper scale to display */
  int k;
  int x_level = (logger->imag_center / logger->y_scale + 1) * height * 0.5;
  int sep_exp = ((int) log10 (logger->x_scale));
  double sep = pow (10, sep_exp) / 3.0;

  double center = logger->real_center / pow (10, sep_exp + 1);
  center = ((int) center) * pow (10, sep_exp + 1);

  char str_buf[1023];

  for (k = -100 ; k <= 100; k++)
  {
    gdouble tick_pos = mps_iteration_logger_x_points_to_coords (logger, center + k * sep);

    if (center + k * sep == 0.0)
      continue;

    cairo_move_to (cr, tick_pos, x_level - 6);
    cairo_line_to (cr, tick_pos, x_level + 6);
    cairo_stroke (cr);

    if (k % 2 == 0)
      {
        sprintf (str_buf, "%1.2e", center + k * sep);
        cairo_text_extents (cr, str_buf, &extents);
        cairo_move_to (cr, tick_pos - extents.width / 2, x_level - extents.height - 2);
        cairo_show_text (cr, str_buf);
      }
  }
}

static void
mps_iteration_logger_draw_y_ticks (MpsIterationLogger * logger, cairo_t * cr)
{
  int width = gtk_widget_get_allocated_width (logger->drawing_area);
  cairo_text_extents_t extents;

  cairo_set_source_rgb (cr, 0.4, 0.4, 0.4);

  /* Get a proper scale to display */
  int k;
  int y_level = (-logger->real_center / logger->x_scale + 1) * width * 0.5;
  int sep_exp = ((int) log10 (logger->y_scale));
  double sep = pow (10, sep_exp) / 3.0;

  double center = logger->imag_center / pow (10, sep_exp + 1);
  center = ((int) center) * pow (10, sep_exp + 1);

  char str_buf[1023];

  for (k = -100 ; k <= 100; k++)
  {
    gdouble tick_pos = mps_iteration_logger_y_points_to_coords (logger, center + k * sep);

    if (center + k * sep == 0.0)
      continue;

    cairo_move_to (cr, y_level - 6, tick_pos);
    cairo_line_to (cr, y_level + 6, tick_pos);
    cairo_stroke (cr);

    if (k % 2 == 0)
      {
        sprintf (str_buf, "%1.2e", center + k * sep);
        cairo_text_extents (cr, str_buf, &extents);
        cairo_move_to (cr, y_level + 6, 
          tick_pos + extents.height / 2);
        cairo_show_text (cr, str_buf);
      }
  }
}

static void
mps_iteration_logger_on_drawing_area_draw (GtkWidget * widget,
                                           cairo_t * cr, MpsIterationLogger * logger)
{
  int width, height;
  double r, g, b;

  if (logger->exit)
  {
    logger->drawing = FALSE;
    return;
  }

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

  /* Draw some indexes on the axes */
  mps_iteration_logger_draw_x_ticks (logger, cr);
  mps_iteration_logger_draw_y_ticks (logger, cr);

  int degree = logger->ctx ? logger->ctx->n : logger->degree;
  mps_approximation ** approximations = logger->ctx ? logger->ctx->root : logger->approximations;

  /* Draw points if present */
  if (approximations)
    {
      int i;
      double x, y;
      cairo_set_source_rgb (cr, 0.9, 0.1, 0.1);

      for (i = 0; i < degree; i++)
        {
          switch (logger->ctx ? logger->ctx->lastphase : mp_phase)
          {
            case mp_phase:
              mpc_get_cplx (approximations[i]->fvalue, approximations[i]->mvalue);

            case dpe_phase:
              cdpe_get_x (approximations[i]->fvalue, approximations[i]->dvalue);

            default:
              x = mps_iteration_logger_x_points_to_coords (logger, 
                cplx_Re (approximations[i]->fvalue));
              y = mps_iteration_logger_y_points_to_coords (logger, 
                cplx_Im (approximations[i]->fvalue));
              break;
          }

          if (approximations[i]->approximated)
          {
            r = 0.1 ; g = 0.9 ; b = 0.1;
          }
          else if (!approximations[i]->again)
          {
            r = 0.1; g = 0.1; b = 0.9;
          }
          else 
          {
            r = 0.9; g = 0.1; b = 0.1;
          }

          /* Check if the user has zommed enough to see the radii */
          if (approximations[i]->frad > 1.3 * MAX (logger->x_scale, logger->y_scale) && false)
            {
              cairo_save (cr);

              cairo_scale (cr, 1.0 / logger->x_scale, 
                1.0 / logger->y_scale);

              cairo_set_source_rgba (cr, r, g, b, 0.3);
              cairo_arc (cr, x, y, approximations[i]->frad, 0, 6.29);
              cairo_fill_preserve (cr);

              cairo_set_source_rgba (cr, r, g, b, 1.0);
              cairo_stroke (cr);

              cairo_restore (cr);
            }
          else
            { 
              cairo_set_source_rgb (cr, r, g, b);
              cairo_arc (cr, x, y, 1.3, 0, 6.29);
              cairo_fill (cr);
            }
        }

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

  logger->drawing = FALSE;
}

static gboolean 
mps_iteration_logger_update_drawing_area (gpointer user_data)
{
  MpsIterationLogger * logger = MPS_ITERATION_LOGGER (user_data);

  /* Fail silently if we are already drawing */
  if (!logger->exit && !logger->drawing)
    {
      logger->drawing = TRUE;
      gtk_widget_queue_draw (logger->drawing_area);
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

  if (!logger->zooming)
    {
      logger->zooming = true;
      logger->zoom_rect_x = event->x;
      logger->zoom_rect_y = event->y;
    }
  else 
    {
      logger->zooming = FALSE;

      logger->real_center = (mps_iteration_logger_x_coords_to_points (logger, event->x) +
        mps_iteration_logger_x_coords_to_points (logger, logger->zoom_rect_x)) / 2;
      logger->imag_center = (mps_iteration_logger_y_coords_to_points (logger, event->y) + 
        mps_iteration_logger_y_coords_to_points (logger, logger->zoom_rect_y)) / 2;

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
  /* We will use this drawing area to draw the computed roots. */
  logger->drawing_area = gtk_drawing_area_new ();

  /* Setup the drawing area to to handle the appropriate events */
  gtk_widget_add_events (logger->drawing_area, GDK_BUTTON_PRESS_MASK);
  g_signal_connect (logger->drawing_area, "button-press-event", 
    G_CALLBACK (mps_iteration_logger_on_da_button_press), logger);

#if GTK_MAJOR_VERSION < 3
  GtkWidget *box = gtk_vbox_new (FALSE, 6);
#else
  GtkWidget *box = gtk_box_new (GTK_ORIENTATION_VERTICAL, 6);
#endif

  gtk_box_pack_start (GTK_BOX (box), logger->drawing_area, TRUE, TRUE, 0);

  /* Create the range adjuster */
  GtkWidget * zoom_out = gtk_button_new_with_label ("Zoom out");
  GtkWidget * zoom_in  = gtk_button_new_with_label ("Reset zoom");

  g_signal_connect (zoom_out, "clicked", G_CALLBACK (mps_iteration_logger_on_zoom_out_clicked),
    logger);
  g_signal_connect (zoom_in, "clicked", G_CALLBACK (mps_iteration_logger_on_reset_zoom_clicked),
    logger);

#if GTK_MAJOR_VERSION < 3
  GtkWidget * zoom_box = gtk_hbox_new (FALSE, 6);
#else
  GtkWidget * zoom_box = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 6);
#endif

  gtk_box_pack_start (GTK_BOX (zoom_box), zoom_in, TRUE, TRUE, 0);
  gtk_box_pack_start (GTK_BOX (zoom_box), zoom_out, TRUE, TRUE, 0);

  gtk_box_pack_start (GTK_BOX (box), zoom_box, FALSE, TRUE, 0);

  gtk_container_add (GTK_CONTAINER (logger), box);

#if GTK_MAJOR_VERSION < 3
  g_signal_connect (logger->drawing_area, "expose-event", 
    G_CALLBACK (mps_iteration_logger_on_expose_event), logger);
#else
  g_signal_connect (logger->drawing_area, "draw", 
    G_CALLBACK (mps_iteration_logger_on_drawing_area_draw), logger);
#endif

  gtk_window_set_title (GTK_WINDOW (logger), "MPSolve iterations");

  /* Start the loop */
  logger->timeout_source = g_timeout_add (1000.0 / 25, 
    mps_iteration_logger_update_drawing_area, logger);

  gtk_widget_set_size_request (GTK_WIDGET (logger), 640, 640);
}

#endif