# coding: utf-8

plt.figure(figsize=[8, 8])
plt.plot(r_0[:, 0], r_0[:, 1], 'bo')

if any(v_0[:, 0]) or any(v_0[:, 1]):
    plt.quiver(r_0[:, 0], r_0[:, 1], v_0[:, 0], v_0[:, 1], scale=velocity_scale,
               pivot='middle', headlength=2, headaxislength=1)
if any(F_0[:, 0]) or any(F_0[:, 1]):
    plt.quiver(r_0[:, 0], r_0[:, 1], F_0[:, 0], F_0[:, 1], scale=force_scale,
               color='red', pivot='middle', headlength=2, headaxislength=1)

plt.xlim(-L, 2 * L)
plt.ylim(-L, 2 * L)

plt.hlines([0, L], -L, 2 * L, linestyles='dashed')
plt.vlines([0, L], -L, 2 * L, linestyles='dashed')

for m in [-L, 0, L]:
    for n in [-L, 0, L]:
        if n is not 0 or m is not 0:
            plt.plot(r_0[:, 0] + n, r_0[:, 1] + m, 'bo', alpha=0.5)
            if any(v_0[:, 0]) or any(v_0[:, 1]):
                plt.quiver(r_0[:, 0] + n, r_0[:, 1] + m, v_0[:, 0], v_0[:, 1],
                           scale=velocity_scale, pivot='middle', alpha=0.5,
                           headlength=2, headaxislength=1)

circle = plt.Circle(r_0[0], cutoff, color='k', fill=False, linestyle='dashdot')
plt.gca().add_artist(circle)

plt.savefig('init.pdf')

# plt.figure(figsize=[8, 8])
# plt.plot(r_0[:, 0], r_0[:, 1], 'bo')
# plt.plot(r_1[:, 0], r_1[:, 1], 'bo')
# if any(v_1[:, 0]) or any(v_1[:, 1]):
#     plt.quiver(r_1[:, 0], r_0[:, 1], v_1[:, 0], v_1[:, 1], scale=velocity_scale,
#                pivot='middle', headlength=2, headaxislength=1)
# if any(F_1[:, 0]) or any(F_1[:, 1]):
#     plt.quiver(r_1[:, 0], r_1[:, 1], F_1[:, 0], F_1[:, 1], scale=force_scale,
#                color='red', pivot='middle', headlength=2, headaxislength=1)

# plt.xlim(0, L)
# plt.ylim(0, L)

# plt.hlines([0, L], -L, 2 * L, linestyles='dashed')
# plt.vlines([0, L], -L, 2 * L, linestyles='dashed')

# for m in [-L, 0, L]:
#     for n in [-L, 0, L]:
#         if n is not 0 or m is not 0:
#             plt.plot(r_1[:, 0] + n, r_1[:, 1] + m, 'bo', alpha=0.5)
#             plt.quiver(r_1[:, 0] + n, r_1[:, 1] + m, v_1[:, 0], v_1[:, 1],
#                        pivot='middle', alpha=0.5, headlength=2, headaxislength=1)

# plt.savefig('first_step.pdf')
