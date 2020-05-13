% Jacobian of implicit kinematic constraints of
% palh1m1IC
% with respect to passive joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% Phi_p [(no of constraints)x(no of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_p = palh1m1IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:56:38
% EndTime: 2020-04-15 19:56:38
% DurationCPUTime: 0.04s
% Computational Cost: add. (28->11), mult. (18->14), div. (0->0), fcn. (18->14), ass. (0->10)
t28 = qJ(8) + qJ(9);
t29 = sin(t28) * pkin(12);
t27 = -qJ(7) + pkin(19);
t25 = pkin(20) + qJ(7) + qJ(2);
t24 = qJ(3) + qJ(4) + pkin(18);
t23 = -qJ(10) + t27;
t22 = cos(t28) * pkin(12);
t21 = cos(t23) * pkin(8);
t20 = sin(t23) * pkin(8);
t1 = [0, -sin(qJ(6)) * pkin(7), pkin(3) * cos(t25), 0, 0, 0, 0, 0, 0; 0, cos(qJ(6)) * pkin(7), pkin(3) * sin(t25), 0, 0, 0, 0, 0, 0; 0, 0, 0, -t29 + sin(qJ(8)) * pkin(2), -t29, 0, 0, 0, 0; 0, 0, 0, t22 - cos(qJ(8)) * pkin(2), t22, 0, 0, 0, 0; pkin(10) * cos(t24), 0, t20 - pkin(4) * sin(t27), 0, 0, t20, 0, 0, 0; pkin(10) * sin(t24), 0, t21 - pkin(4) * cos(t27), 0, 0, t21, 0, 0, 0; 0, 1, -1, 0, 0, 0, -1, 0, 0; 0, 0, 0, -1, -1, 0, 0, 1, 0; -1, 0, 1, 0, 0, 1, 0, 0, -1;];
Phi_p = t1;
