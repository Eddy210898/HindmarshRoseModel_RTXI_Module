/*
 * Copyright (C) 2011 Georgia Institute of Technology, University of Utah,
 * Weill Cornell Medical College
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This is a template implementation file for a user module derived from
 * DefaultGUIModel with a custom GUI.
 */

#include "HindmarshRoseModelPlugin.h"
#include <iostream>
#include <main_window.h>

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
  return new HindmarshRoseModelPlugin();
}

static DefaultGUIModel::variable_t vars[] = {
    {
        "I",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },

    {
        "a",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "b",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "c",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "d",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "r",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "s",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },

    {
        "x0",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "y0",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "z0",
        "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },

    {
        "Vm (mV)",
        "Membrane potential (in mV)",
        DefaultGUIModel::OUTPUT,
    },

    {
        "Isyn",
        "",
        DefaultGUIModel::INPUT,
    },
    {
        "v",
        "",
        DefaultGUIModel::STATE,
    },
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

HindmarshRoseModelPlugin::HindmarshRoseModelPlugin(void)
    : DefaultGUIModel("HindmarshRoseModelPlugin with Custom GUI", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>HindmarshRoseModelPlugin:</b><br>QWhatsThis description.</p>");
  DefaultGUIModel::createGUI(vars,
                             num_vars); // this is required to create the GUI
  customizeGUI();
  initParameters();
  update(INIT); // this is optional, you may place initialization code directly
                // into the constructor
  refresh();    // this is required to update the GUI with parameter and state
                // values
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

HindmarshRoseModelPlugin::~HindmarshRoseModelPlugin(void)
{
}
double HindmarshRoseModelPlugin::getXValue(int kValue, double x, double dX)
{
  if (kValue == 1)
  {
    return x;
  }
  else if (kValue == 2 || kValue == 3)
  {
    return x + (0.5 * dX);
  }
  else if (kValue == 4)
  {
    return x + dX;
  }
  else
  {
    throw "K value not recognize";
  }
}

double HindmarshRoseModelPlugin::getYValue(int kValue, double y, double dX, double k /*= 0*/)
{
  if (kValue == 1)
  {
    return y;
  }
  else if (kValue == 2 || kValue == 3)
  {
    return y + (0.5 * k * dX);
  }
  else if (kValue == 4)
  {
    return y + (k * dX);
  }
  else
  {
    throw "K value not recognize";
  }
}

double HindmarshRoseModelPlugin::getNextRungeKuta(double Xo, double Yo, double dX, double F(double, double, double[]), double args[])
{
  double x = getXValue(1, Xo, dX);
  double y = getYValue(1, Yo, dX);
  double kA = F(x, y, args);
  x = getXValue(2, Xo, dX);
  y = getYValue(2, Yo, dX, kA);
  double kB = F(x, y, args);
  x = getXValue(3, Xo, dX);
  y = getYValue(3, Yo, dX, kA);
  double kC = F(x, y, args);
  x = getXValue(4, Xo, dX);
  y = getYValue(4, Yo, dX, kA);
  double kD = F(x, y, args);
  double kT = kA + kB + kC + kD;
  double newValToAdd = (dX * kT) / 6;
  return Yo + newValToAdd;
}
void HindmarshRoseModelPlugin::hindmarshRoseStep(double xO, double yO, double zO, double xOr, double t, double dT, double I, double a, double b, double c, double d, double r, double s)
{
  double arg[] = {
      a,   //0
      b,   //1
      c,   //2
      d,   //3
      xO,  //4
      xOr, //5
      yO,  //6
      zO,  //7
      I,   //8
      r,   //9
      s    //10
  };
  double xF = getNextRungeKuta(
      t, xO, dT, [](double t, double x, double args[]) {
        double sE = args[0] * x * x * x;
        double tE = args[1] * x * x;
        return args[6] - sE + tE + args[8] - args[7];
      },
      arg);

  double yF = getNextRungeKuta(
      t, yO, dT, [](double t, double y, double args[]) {
        double sE = args[3] * args[4] * args[4];
        return args[2] - sE - y;
      },
      arg);
  double zF = getNextRungeKuta(
      t, zO, dT, [](double t, double z, double args[]) {
        double sE = args[4] - args[5];
        double tE = args[10] * sE;
        double fE = tE - args[7];
        return args[9] * fE;
      },
      arg);

  x = xF;
  y = yF;
  z = zF;
}

void HindmarshRoseModelPlugin::execute(void)
{
  hindmarshRoseStep(x, y, z, xO, period, dt, I, a, b, c, d, r, s);
  output(0) = x;
  return;
}

void HindmarshRoseModelPlugin::initParameters(void)
{
  a = 1;
  b = 3;
  c = 1;
  d = 5;
  r = 0.001;
  s = 1;
  I = 1;
  x = -1.5;
  xO = x;
  y = 1 - (5 * x * x);
  z = s * (x + 1.6);
}

void HindmarshRoseModelPlugin::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag)
  {
  case INIT:
    period = RT::System::getInstance()->getPeriod() * 1e-8; // ms

    setParameter("I", I);

    setParameter("a", a);
    setParameter("b", b);
    setParameter("c", c);
    setParameter("d", d);
    setParameter("r", r);
    setParameter("s", s);

    setParameter("x0", x);
    setParameter("y0", y);
    setParameter("z0", z);

    setState("v", x);
    break;

  case MODIFY:
    a = getParameter("a").toDouble();
    b = getParameter("b").toDouble();
    c = getParameter("c").toDouble();
    d = getParameter("d").toDouble();
    r = getParameter("r").toDouble();
    s = getParameter("s").toDouble();

    I = getParameter("I").toDouble();

    x = getParameter("x0").toDouble();
    y = getParameter("y0").toDouble();
    z = getParameter("z0").toDouble();
    break;

  case UNPAUSE:
    break;

  case PAUSE:
    break;

  case PERIOD:
    period = RT::System::getInstance()->getPeriod() * 1e-8; // ms
    break;

  default:
    break;
  }
}

void HindmarshRoseModelPlugin::customizeGUI(void)
{
  /*QGridLayout *customlayout = DefaultGUIModel::getLayout();

  QGroupBox *button_group = new QGroupBox;

  QPushButton *abutton = new QPushButton("Button A");
  QPushButton *bbutton = new QPushButton("Button B");
  QHBoxLayout *button_layout = new QHBoxLayout;
  button_group->setLayout(button_layout);
  button_layout->addWidget(abutton);
  button_layout->addWidget(bbutton);
  QObject::connect(abutton, SIGNAL(clicked()), this, SLOT(aBttn_event()));
  QObject::connect(bbutton, SIGNAL(clicked()), this, SLOT(bBttn_event()));

  customlayout->addWidget(button_group, 0, 0);
  setLayout(customlayout);*/
}
