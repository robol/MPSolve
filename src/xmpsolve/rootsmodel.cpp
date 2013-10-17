#include "rootsmodel.h"
#include <QDebug>

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

    role_names.insert(RADIUS, "radius");
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
        rdpe_t root_module;
        mpc_rmod (root_module, m_roots.at(i)->value);
        int digits = (rdpe_Esp (root_module) - rdpe_Esp (m_roots.at(i)->radius)) / LOG2_10 + 1;
        char * buffer = NULL;

        switch (role)
        {
            case SHORT_APPROXIMATION:
                digits = 4;
                // fallthrough
            case Qt::DisplayRole:
	        buffer = new char[2 * (digits + 15)];
                if (m_roots[i]->get_imag_part() > 0)
                    gmp_sprintf (buffer, "%.*Fe + %.*Fei", digits, mpc_Re (m_roots[i]->value),
                                 digits, mpc_Im (m_roots[i]->value));
                else
                    gmp_sprintf (buffer, "%.*Fe %.*Fei", digits, mpc_Re (m_roots[i]->value),
                                 digits, mpc_Im (m_roots[i]->value));
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
