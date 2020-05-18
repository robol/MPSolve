#include "rootsmodel.h"
#include <iostream>
#include <QDebug>

static double LOG2_10 = log(10) / log(2);

namespace xmpsolve {

RootsModel::RootsModel(QObject *parent) :
    QAbstractListModel(parent)
{
    m_length = 0;
    m_marked_root = -1;
}

QHash<int, QByteArray>
RootsModel::roleNames() const
{
    QHash<int, QByteArray> role_names;

    role_names.insert(RADIUS, "inclusion_radius");
    role_names.insert(STATUS, "status");
    role_names.insert(SHORT_APPROXIMATION, "short_approximation");
    role_names.insert(ROOT, "root");
    role_names.insert(MARKED, "marked");

    return role_names;
}

int
RootsModel::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return m_length;
}

QVariant
RootsModel::data(const QModelIndex &index, int role) const
{
    int i = index.row();

    if (role != Qt::DisplayRole && role < Qt::UserRole)
        return QVariant();

    if (i < 0 || i > m_length)
        return QVariant();
    else
    {
	cdpe_t croot;
	mpcf_get_cdpe (croot, m_roots.at(i)->value);
	int real_digits = (rdpe_Esp (cdpe_Re (croot)) - rdpe_Esp (m_roots.at(i)->radius)) / LOG2_10 - 1;
	int imag_digits = (rdpe_Esp (cdpe_Im (croot)) - rdpe_Esp (m_roots.at(i)->radius)) / LOG2_10 - 1;

        char * buffer = NULL;

        int imag_digits_size = (imag_digits < 0) ? 2 : imag_digits;
        int real_digits_size = (real_digits < 0) ? 2 : real_digits;

        switch (role)
        {
            case SHORT_APPROXIMATION:
                real_digits = (real_digits < 4) ? real_digits : 4;
                imag_digits = (imag_digits < 4) ? imag_digits : 4;
	        // fallthrough
            case Qt::DisplayRole:
	        buffer = new char[real_digits_size + imag_digits_size + 20];

		if (imag_digits <= 0 && real_digits <= 0)
		  gmp_sprintf (buffer, "0.0");
		else if (imag_digits <= 0)
		  gmp_sprintf (buffer, "%.*Fe", real_digits, mpcf_Re (m_roots[i]->value));
		else if (real_digits <= 0)
		  gmp_sprintf (buffer, "%.*Fei", imag_digits, mpcf_Im (m_roots[i]->value));
		else 
		  {
		    if (m_roots[i]->get_imag_part() > 0)
		      gmp_sprintf (buffer, "%.*Fe + %.*Fei", real_digits, mpcf_Re (m_roots[i]->value),
				   imag_digits, mpcf_Im (m_roots[i]->value));
		    else
		      gmp_sprintf (buffer, "%.*Fe %.*Fei", real_digits, mpcf_Re (m_roots[i]->value),
				   imag_digits, mpcf_Im (m_roots[i]->value));
		  }

                return QString(buffer);
                break;

            case STATUS:
                return QString(MPS_ROOT_STATUS_TO_STRING (m_roots[i]->status));

            case RADIUS:
                return QString("%1").arg(m_roots[i]->get_radius());

            case ROOT:
                return QVariant::fromValue((void*) m_roots[i]);

            case MARKED:
                return i == m_marked_root;

            default:
                qDebug() << "Invalid role";
                return QVariant();
        }
    }
}

void
RootsModel::setRoots(QList<Root *> roots)
{
    beginResetModel();

    m_length = 0;
    m_roots = roots;
    m_length = roots.length();

    endResetModel();
}

void
RootsModel::markRoot(int i)
{
    int oldMarkedRoot = m_marked_root;

    if (i >= -1 && i < m_length)
        m_marked_root = i;
    else
        m_marked_root = -1;

    if (oldMarkedRoot != -1)
        dataChanged(index(oldMarkedRoot), index(oldMarkedRoot));
    if (m_marked_root != -1)
        dataChanged(index(m_marked_root), index(m_marked_root));
}



} // namespace xmpsolve
