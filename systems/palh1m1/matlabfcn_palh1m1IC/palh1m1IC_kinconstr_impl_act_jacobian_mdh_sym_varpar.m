% Jacobian of implicit kinematic constraints of
% palh1m1IC
% with respect to active joint coordinates
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
% Phi_a [(no of constraints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_a = palh1m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:56:38
% EndTime: 2020-04-15 19:56:38
% DurationCPUTime: 0.02s
% Computational Cost: add. (14->7), mult. (10->10), div. (0->0), fcn. (10->10), ass. (0->4)
t9 = qJ(3) + pkin(17);
t8 = pkin(20) + qJ(7) + qJ(2);
t7 = qJ(3) + qJ(4) + pkin(18);
t1 = [0, pkin(3) * cos(t8) + cos(qJ(2)) * pkin(1), 0, 0; 0, pkin(3) * sin(t8) + sin(qJ(2)) * pkin(1), 0, 0; 0, 0, pkin(6) * cos(t9), 0; 0, 0, pkin(6) * sin(t9), 0; 0, 0, pkin(10) * cos(t7) + cos(qJ(3)) * pkin(5), 0; 0, 0, pkin(10) * sin(t7) + sin(qJ(3)) * pkin(5), 0; 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, -1, 0;];
Phi_a = t1;
