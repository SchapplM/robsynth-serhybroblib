% Implicit kinematic constraints of
% palh1m1IC
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
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = palh1m1IC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_kinconstr_impl_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_kinconstr_impl_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:56:38
% EndTime: 2020-04-15 19:56:38
% DurationCPUTime: 0.06s
% Computational Cost: add. (53->33), mult. (23->23), div. (0->0), fcn. (20->20), ass. (0->7)
t6 = qJ(8) + qJ(9);
t5 = qJ(3) + pkin(17);
t4 = -qJ(7) + pkin(19);
t3 = pkin(20) + qJ(7) + qJ(2);
t2 = qJ(3) + qJ(4) + pkin(18);
t1 = -qJ(10) + t4;
t7 = [cos(qJ(6)) * pkin(7) - pkin(14) + pkin(3) * sin(t3) + sin(qJ(2)) * pkin(1) - pkin(15); sin(qJ(6)) * pkin(7) - pkin(16) - pkin(3) * cos(t3) - cos(qJ(2)) * pkin(1); pkin(6) * sin(t5) + pkin(1) + cos(t6) * pkin(12) - cos(qJ(8)) * pkin(2); -pkin(6) * cos(t5) + sin(t6) * pkin(12) - sin(qJ(8)) * pkin(2); pkin(10) * sin(t2) + sin(qJ(3)) * pkin(5) + cos(t1) * pkin(8) - pkin(4) * cos(t4); -pkin(10) * cos(t2) - cos(qJ(3)) * pkin(5) - sin(t1) * pkin(8) + pkin(4) * sin(t4); qJ(6) - qJ(11) + pi / 0.2e1 - t3; qJ(12) + 0.3e1 / 0.2e1 * pi + t5 - t6; -qJ(13) + 0.5e1 / 0.2e1 * pi - t1 - t2;];
h = t7;
